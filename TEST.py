from A_modules.shared.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory


from A_modules.shared.utils.utils import check_directory_exists

test = ConverterFactory().get_converter("pdb", "mol2")
mol2_file = test.run("input/solvents/pdb/hexane.pdb", "TEST", verbose=True)
import pandas as pd

print("!!!" + mol2_file)
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)

parameterizer = ACPYPEParameterizer()
file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)

parametized = parameterizer.run(mol2_file, "TEST", file_config, verbose=True)

from A_modules.shared.pdb_validation.gromacs_validator import GROMACSPDBValidator

updated_path = "input/solvents/pdb/hexane.pdb"

# gromacsvalidator = GROMACSPDBValidator()
# updated_path = "TEST/edited_pdb.pdb"
# gromacsvalidator.validate(
#    "input/solvents/pdb/hexane.pdb",
#    output_file_path=updated_path,
#    padding=0.5 * 10**-10,
# )
gro = parametized.gro_path
top = parametized.top_path

import numpy as np
from typing import List


def calculate_minimum_box_size(
    atom_data: List[List[float]], padding: float = 0.1
) -> List[float]:
    """
    Helper function to calculate the minimum bounding box dimensions for the molecule
    based on the atom coordinates.

    Args:
        atom_data (List[List[float]]): List of parsed atom data, where each entry contains [x, y, z] coordinates.
        padding (float): Padding to add to each dimension (in nm).

    Returns:
        List[float]: [x, y, z] dimensions for the minimum bounding box in nm.
    """
    if not atom_data:
        raise ValueError("Atom data is empty. Cannot calculate bounding box size.")

    # Extract coordinates
    coordinates = np.array([[atom[4], atom[5], atom[6]] for atom in atom_data])

    # Calculate min and max for each dimension
    min_coords = coordinates.min(axis=0)
    max_coords = coordinates.max(axis=0)

    # Add padding and calculate box size
    box_size = (max_coords - min_coords) + 2 * padding
    return box_size.tolist()


def calculate_minimum_box_size_from_df(
    atom_df: pd.DataFrame, padding: float = 0.1
) -> List[float]:
    """
    Helper function to calculate the minimum bounding box dimensions for the molecule
    based on the atom coordinates in a DataFrame.

    Args:
        atom_df (pd.DataFrame): DataFrame containing atom data, including 'X', 'Y', and 'Z' columns for coordinates.
        padding (float): Padding to add to each dimension (in nm).

    Returns:
        List[float]: [x, y, z] dimensions for the minimum bounding box in nm.
    """
    if atom_df.empty:
        raise ValueError("Atom data is empty. Cannot calculate bounding box size.")

    # Ensure the required columns are present
    required_columns = ["X", "Y", "Z"]
    if not all(col in atom_df.columns for col in required_columns):
        raise ValueError(
            f"DataFrame must contain columns: {', '.join(required_columns)}"
        )

    # Calculate min and max for each dimension
    min_coords = atom_df[["X", "Y", "Z"]].min().values
    max_coords = atom_df[["X", "Y", "Z"]].max().values

    # Add padding and calculate box size
    box_size = (max_coords - min_coords) + 2 * padding
    return box_size.tolist()


from A_modules.atomistic.gromacs.parser.gromacs_parser import (
    GromacsParser,
)
from A_modules.atomistic.gromacs.parser.registries.handler_registry import (
    handler_registry,
)


file_splitter = GromacsParser()
sections = file_splitter.parse(gro)
test_section = next(iter(sections.values()))
tester = handler_registry.get_handler(test_section.handler_name)()
tester.process(test_section)

print(tester.content)
box_size = calculate_minimum_box_size_from_df(tester.content, padding=0.1)
print(f"Calculated box size: {box_size} nm")
from A_modules.atomistic.gromacs.commands.editconf import Editconf

# editconf = Editconf()
# output_box_gro_path = editconf.run(
#    gro,
#    "TEST",
#    [5, 5, 5],
#    output_name="edited_box.gro",
# )

# from A_modules.atomistic.gromacs.commands.solvate import Solvate

# print(type(updated_path), type(output_box_gro_path), type(top))
# print(output_box_gro_path)
# solvate = Solvate()
# solvate.run(
#    output_box_gro_path,
#    gro,
#    top,
#    "TEST",
#    output_name="solvated.gro",
# )

# from gromacs.solvation.NEW_solvent_insertion import SolventInsertion

# from A_modules.atomistic.gromacs.commands.solvent_insertion import SolventInsertion

# solvent_inserter = SolventInsertion()
# solvent_inserter.run(
#    updated_path,
#    "TEST",
#    660,
#    85,
#    tolerance=0.1,
#    box_size_nm=[5, 5, 5],
##)


# NEED TO RUN SOLVATE FIRST, AND THEN DEAL WITH INSERTION
