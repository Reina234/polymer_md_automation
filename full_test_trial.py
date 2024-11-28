from preprocessing.parsers.pdb_parser import PDBParser
from preprocessing.calculation_utils import (
    calculate_minimum_box_size,
    calculate_density,
    calculate_volume_for_desired_density,
    scale_box_to_desired_volume,
    calculate_num_particles,
)
from preprocessing.validators.pdb_validators.gromacs_pdb_validator import (
    GROMACSPDBValidator,
)
from data_models.solvent import Solvent
from preprocessing.metadata_tracker import MetadataTracker
from gromacs.solvation.polymer_resize_box import PolmerBoxResize
from gromacs.solvation.solvent_box_creator import SolventInsertion
from config.constants import LengthUnits
from gromacs.solvation.solvate import Solvate
from gromacs.gromacs_utils import prepare_topol_file, reformat_topol_file
from gromacs.solvation.ion_adder import IonAdder
from preprocessing.file_converters.obabel_converter import OpenBabelConverter
from preprocessing.metadata_tracker import MetadataTracker
from preprocessing.parameterizers.acpype_parameterizer import ACPYPEParameterizer
from preprocessing.validators.pdb_validators.base_pdb_validator import BasePDBValidator
from preprocessing.parsers.pdb_parser import PDBParser
from preprocessing.parsers.itp_parser import ITPParser
from config.paths import SOLVENT_ITP_DIR, FilesToExtract
from preprocessing.solvent_atomtypes_manager import AtomtypesManager
from gromacs.gromacs_utils import add_atomtypes_to_topology
from gromacs.equilibriation.energy_minimizer import EnergyMinimizer

solvent_files = [FilesToExtract.ITP]

metadata_tracker = MetadataTracker()
pdbtomol2converter = OpenBabelConverter(metadata_tracker)
solvent_pdb_path = "input/solvents/pdb/hexane.pdb"
solvent = Solvent("Hexane", 86.18, 660, solvent_pdb_path, "TMZK")
solvent_mol2_path = pdbtomol2converter.convert(solvent.pdb_path, "temp")

solvent_acpype = ACPYPEParameterizer(metadata_tracker, solvent.pdb_molecule_name)
solvent_itp = solvent_acpype.parameterize(
    solvent_mol2_path,
    SOLVENT_ITP_DIR,
    solvent_files,
    "hexane",
)

solvent_itp = solvent_itp[0]
manager = AtomtypesManager()

manager.extract_and_store_atomtypes(solvent_itp, "hexane")
solvent_itp_content = ITPParser.read_file(solvent_itp)
no_atomtypes = ITPParser.remove_section(solvent_itp_content, "atomtypes")
ITPParser.save_file(solvent_itp, no_atomtypes)

polymer_pdb = "styrene.pdb"
polymer_mol2_path = pdbtomol2converter.convert(polymer_pdb, "output/TRIAL")
acpype = ACPYPEParameterizer(metadata_tracker)

polymer_files = [FilesToExtract.GRO, FilesToExtract.ITP, FilesToExtract.TOP]
acpype_path = "output/TRIAL/acpype_output"
files = acpype.parameterize(polymer_mol2_path, acpype_path, polymer_files, posre=True)
polymer_gro = files[0]
polymer_itp = files[1]
polymer_top = files[2]


topol_file = prepare_topol_file("output/TRIAL/acpype_output/POLY_GMX.top", "TRIAL")

forcefield = "amber99sb-ildn.ff/forcefield.itp"
topol_file = reformat_topol_file(
    topol_file, polymer_itp, solvent_itp, forcefield, "output/TRIAL/gromacs"
)

add_atomtypes_to_topology("hexane", polymer_itp)

gromacsvalidator = GROMACSPDBValidator(metadata_tracker)
updated_solvent_pdb = "hexane_edited.pdb"
gromacsvalidator.validate(solvent.pdb_path, output_file_path=updated_solvent_pdb)
box_creator = PolmerBoxResize(metadata_tracker)
box_file = box_creator.run(polymer_gro, "TRIAL")

solvent_box = SolventInsertion(metadata_tracker)
solvent_box_file = solvent_box.run(
    updated_solvent_pdb, "TRIAL", solvent.density, solvent.molecular_weight
)

print(solvent_box_file, box_file, topol_file)
solute_box = Solvate(metadata_tracker)
solvated_box = solute_box.run(solvent_box_file, box_file, topol_file, "TRIAL")

ion_adder = IonAdder(metadata_tracker)
added_ion_box = ion_adder.run(solvated_box, topol_file, "TRIAL")

energy_minimizer = EnergyMinimizer(metadata_tracker)
minimized_box = energy_minimizer.run(added_ion_box, topol_file, "TRIAL")

##############
# consider separating into solvation and equilibriation
