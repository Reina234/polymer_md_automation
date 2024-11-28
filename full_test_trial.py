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

metadata_tracker = MetadataTracker()
pdbtomol2converter = OpenBabelConverter(metadata_tracker)
solvent_pdb_path = "input/solvents/pdb/hexane.pdb"
solvent = Solvent("Hexane", 86.18, 660, solvent_pdb_path, "TMZK")
solvent_mol2_path = pdbtomol2converter.convert(solvent.pdb_path, "temp")

solvent_acpype = ACPYPEParameterizer(metadata_tracker, solvent.pdb_molecule_name)
solvent_itp = solvent_acpype.parameterize(
    solvent_mol2_path, solvent.pdb_molecule_name, "gromacs", is_solvent=True
)

trial_output_name = "trial_output"

polymer_mol2_path = pdbtomol2converter.convert("styrene.pdb")

acpype = ACPYPEParameterizer(metadata_tracker)
dir = acpype.parameterize(polymer_mol2_path, trial_output_name)
forcefield = "amber99sb-ildn.ff/forcefield.itp"

topol_file = prepare_topol_file(
    "output/trial_output/acpype_output/POLY_GMX.top", "trial_output"
)
topol_file = reformat_topol_file(
    topol_file,
    "output/test_new_struct/acpype_output/POLY_GMX.itp",
    solvent_itp,
    forcefield,
)
