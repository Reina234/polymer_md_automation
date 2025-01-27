from modules.rdkit.polymer_builders.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 5, "rdkit_test2", overwrite=False, save=False)
# print(generator.map)


def itp_to_dataframe(itp_file: str):
    parser = GromacsParser()
    sections = parser.parse("rdkit_test/POLY_GMX.itp")
    atom_section = sections["data_atoms"]
    data_handler = DataHandler()
    # data_handler.section = atom_section
    data_handler.process(atom_section)
    return data_handler.content


def verify_and_generate_mapping(
    mapping_data: List[dict], itp_df: pd.DataFrame, output_file: str
):
    """
    Verify and generate the PyCGTOOL mapping file, ensuring alignment with .itp headers.

    :param mapping_data: List of mapping dictionaries from the polymer generator.
    :param itp_df: Parsed .itp file as a DataFrame with expected headers.
    :param output_file: Path to the mapping file.
    """
    # Ensure .itp DataFrame contains expected columns
    required_columns = {"nr", "type", "resi", "res", "atom"}
    if not required_columns.issubset(itp_df.columns):
        raise ValueError(
            f".itp file is missing required columns. Expected: {required_columns}, "
            f"but got: {set(itp_df.columns)}"
        )

    # Ensure atom_index is integer for proper merging
    itp_df["nr"] = itp_df["nr"].astype(int)

    # Create a lookup dictionary for atom names based on `nr` (atom index)
    atom_name_lookup = itp_df.set_index("nr")["atom"].to_dict()

    # Generate mapping file
    with open(output_file, "w") as f:
        f.write("[MAPPING]\n")
        for entry in mapping_data:
            atom_index = entry["atom_index"]
            bead_type = entry["bead_type"]
            unique_name = entry["unique_name"]

            # Lookup atom name from .itp file
            atom_name = atom_name_lookup.get(atom_index, "UNKNOWN")

            # Write mapping entry
            f.write(f"{unique_name} {bead_type} {atom_name}\n")

    print(f"Mapping file '{output_file}' generated successfully.")


itp_df = itp_to_dataframe("rdkit_test/POLY_GMX.itp")
verify_and_generate_mapping(generator.pycg_map, itp_df, "mapping2.txt")
print(generator.cg_map)
