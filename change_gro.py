from typing import List, Optional
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.utils.atomistic.file_utils import calculate_minimum_box_size_from_df
from modules.gromacs.parsers.handlers.gro_handler import GroHandler
import re


def add_box_dim(
    gro_file, padding: float = 0.1, parser=GromacsParser(), gro_handler=GroHandler()
) -> List[float]:
    sections = parser.parse(gro_file)
    first_key = next(iter(sections))  # Get the first key
    gro_section = sections[first_key]
    gro_handler.process(gro_section)
    box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding)
    gro_handler.box_dimensions = box_size
    sections[first_key] = gro_handler.export()

    parser.export(sections, gro_file)

    return gro_file


def shift_atom_labels(
    gro_file,
    output_file: Optional[str] = None,
    padding: float = 0.1,
    parser=GromacsParser(),
    gro_handler=GroHandler(),
) -> List[float]:
    sections = parser.parse(gro_file)
    first_key = next(iter(sections))  # Get the first key
    gro_section = sections[first_key]
    gro_handler.process(gro_section)
    shifted_content = shift_atom_labels_df(gro_handler.content)
    gro_handler.content = shifted_content
    sections[first_key] = gro_handler.export()

    if not output_file:
        output_file = gro_file
    parser.export(sections, output_file)

    return output_file


def shift_atom_labels_df(df, column_name="Atom Name"):
    """
    Shifts the indices of atoms in a DataFrame column by decrementing numerical suffixes.

    Args:
        df (pd.DataFrame): The input DataFrame containing the atom names.
        column_name (str): The column in which atom names should be modified.

    Returns:
        pd.DataFrame: The updated DataFrame with shifted atom names.
    """

    def decrement_atom_label(atom_name):
        """Shift numerical suffix down by 1 if present."""
        match = re.match(r"([A-Za-z]+)(\d+)", atom_name)
        if match:
            base, num = match.groups()
            new_num = int(num) - 1
            return (
                f"{base}{new_num}" if new_num >= 0 else base
            )  # Avoid negative indices
        return atom_name  # Return unchanged if no number is found

    df[column_name] = df[column_name].apply(decrement_atom_label)
    print(df)
    return df


file = shift_atom_labels("rdkit_test3/c=cc1ccccc1_10.gro", "manual_gro.gro")
print(file)
box_file = add_box_dim(file)
