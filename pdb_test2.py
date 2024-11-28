from preprocessing.parsers.pdb_parser import PDBParser
from preprocessing.pdb_utils import (
    calculate_minimum_box_size,
    calculate_density,
    calculate_volume_for_desired_density,
    scale_box_to_desired_volume,
)
from preprocessing.validators.pdb_validators.gromacs_pdb_validator import (
    GROMACSPDBValidator,
)
from data_models.solvent import Solvent
from preprocessing.metadata_tracker import MetadataTracker

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

gromacsvalidator.validate("TEST.pdb", solvent, output_file_path="TEST2.pdb")
