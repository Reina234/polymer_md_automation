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
from gromacs.solvation.polymer_box_creator import PolmerBoxResize
from gromacs.solvation.solvent_box_creator import SolventInsertion
from config.constants import LengthUnits
from gromacs.solvation.solvate import Solvate
from gromacs.gromacs_utils import prepare_topol_file

solvent = Solvent(
    "Hexane",
    86.18,
    0.660,
    "input/solvents/itp/hexane.itp",
    "input/solvents/pdb/hexane.pdb",
    "gromos54a7.ff",
)
# Initialize parser
# parser = PDBParser(file_path="styrene.pdb")

# Extract atom coordinates and calculate minimum box size
# atom_coords = parser.get_atom_coordinates()
# min_box_size = calculate_minimum_box_size(atom_coords)

# print(min_box_size)
# Calculate density based on extracted box dimensions
# density = calculate_density(molecular_weight=10, box_dimensions=min_box_size)

# print(f"Minimum box size: {min_box_size}")
# print(f"Calculated density: {density:.6f} g/cm³")


# required_volume = calculate_volume_for_desired_density(
#    molecular_weight=10, desired_density=1.0
# )
# print(required_volume)


# box_dim = scale_box_to_desired_volume(min_box_size, required_volume)
# print(box_dim)

gromacsvalidator = GROMACSPDBValidator(metadata_tracker=MetadataTracker())

gromacsvalidator.validate("TEST.pdb", output_file_path="TEST2.pdb")

solvent = Solvent(
    "Hexane",
    86.18,
    0.660 * 1000,
    "input/solvents/itp/hexane.itp",
    "TEST2.pdb",
    "gromos54a7.ff",
)

box_creator = PolmerBoxResize()
box_file = box_creator.run(
    "output/test_new_struct/acpype_output/POLY_GMX.gro", "test_new_struct"
)
############## need to sort out file locations
solvent_box = SolventInsertion()
solvent_file = solvent_box.run(
    solvent.pdb_path, "test_new_struct", solvent.density, solvent.molecular_weight
)

topol_file = prepare_topol_file(
    "output/test_new_struct/acpype_output/POLY_GMX.top", "test_new_struct"
)

solute_box = Solvate()
solvated_box = solute_box.run(solvent_file, box_file, topol_file, "test_new_struct")
#    "test_new_struct",
# )
