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
from gromacs.solvation.NEW_solvent_insertion import SolventInsertion
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
from gromacs.equilibriation.temperature_equilibriation import TemperatureEquilibration
from gromacs.equilibriation.pressure_equilibriation import PressureEquilibration

solvent_files = [FilesToExtract.ITP, FilesToExtract.TOP]

metadata_tracker = MetadataTracker()
pdbtomol2converter = OpenBabelConverter(metadata_tracker)
solvent_pdb_path = "input/solvents/pdb/hexane.pdb"
solvent = Solvent("Hexane", 86.18, 660, solvent_pdb_path, "TMZK")
solvent_mol2_path = pdbtomol2converter.convert(solvent.pdb_path, "temp")

solvent_acpype = ACPYPEParameterizer(metadata_tracker, solvent.pdb_molecule_name)
solvent_files = solvent_acpype.parameterize(
    solvent_mol2_path,
    "solvent_trial",
    solvent_files,
    "hexane",
)

solvent_itp = solvent_files[0]
solvent_topol = solvent_files[1]
# manager = AtomtypesManager()

# manager.extract_and_store_atomtypes(solvent_itp, "hexane")
# solvent_itp_content = ITPParser.read_file(solvent_itp)
# no_atomtypes = ITPParser.remove_section(solvent_itp_content, "atomtypes")
# ITPParser.save_file(solvent_itp, no_atomtypes)


gromacsvalidator = GROMACSPDBValidator(metadata_tracker)
updated_solvent_pdb = "hexane_edited.pdb"
gromacsvalidator.validate(solvent.pdb_path, output_file_path=updated_solvent_pdb)

solvent_box = SolventInsertion(metadata_tracker)
solvent_box_file = solvent_box.run(
    updated_solvent_pdb,
    "solvent_trial_2",
    solvent.density,
    solvent.molecular_weight,
    box_size_nm=[5, 5, 5],
)
