from modules.shared.utils.dataframe_utils import (
    dataframe_not_empty_check,
)
from modules.shared.utils.file_utils import (
    file_exists_check_wrapper,
    file_type_check_wrapper,
    check_directory_exists,
    prepare_output_file_path,
    add_suffix_to_filename,
    copy_file,
)
from modules.atomistic.gromacs.parser.handlers.includes_handler import IncludesHandler
from collections import OrderedDict
from modules.atomistic.gromacs.commands.insert_molecules import InsertMolecules
from config.paths import TEMP_DIR
from config.constants import MassUnits, LengthUnits
from modules.shared.utils.calculation_utils import calculate_num_particles
import pandas as pd
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler
from typing import List, Optional, Dict, Union
from data_models.solvent import Solvent
from modules.atomistic.gromacs.commands.solvate import Solvate
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
import logging
import os
from modules.atomistic.gromacs.commands.editconf import Editconf
from data_models.output_types import GromacsPaths
from modules.atomistic.gromacs.parser.data_models.section import Section

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def delete_all_include_sections(sections: OrderedDict) -> OrderedDict:
    """
    Deletes all sections with construct_name set to 'include'.

    Args:
        sections (OrderedDict): The parsed sections from the topology file.

    Returns:
        OrderedDict: The updated sections with all 'include' sections removed.
    """
    return OrderedDict(
        (key, section)
        for key, section in sections.items()
        if section.construct_name != "include"
    )


def calculate_molecule_counts(
    gro_handler: GroHandler,
    residue_name_col: str = "Residue Name",
    residue_number_col: str = "Residue Number",
) -> pd.DataFrame:
    """
    Calculates the number of molecules for each unique residue name by counting unique residue numbers,
    preserving the order of first occurrence.

    :param gro_handler: Input GroHandler containing residue information.
    :param residue_name_col: Column name for residue names.
    :param residue_number_col: Column name for residue numbers.
    :return: A dataframe with residue names and their corresponding molecule counts, in the order of first occurrence.
    """
    dataframe = gro_handler.content
    print(gro_handler.content)

    # Ensure residue names appear in their first occurrence order
    residue_name_order = dataframe[residue_name_col].drop_duplicates()
    dataframe[residue_name_col] = pd.Categorical(
        dataframe[residue_name_col], categories=residue_name_order, ordered=True
    )

    # Group by residue name and count unique residue numbers
    molecule_counts = (
        dataframe.groupby(residue_name_col, observed=True)[residue_number_col]
        .nunique()
        .reset_index()
        .rename(
            columns={
                residue_name_col: "Residue Name",
                residue_number_col: "Number of Molecules",
            }
        )
    )

    return molecule_counts


def replace_value_in_dataframe(
    dataframe: pd.DataFrame, target_value: str, replacement_value: str
) -> pd.DataFrame:
    """
    Replace all occurrences of a specific value in the DataFrame with a new value.

    :param dataframe: The input DataFrame where values will be replaced.
    :param target_value: The value to be replaced.
    :param replacement_value: The value to replace with.
    :return: The updated DataFrame with the values replaced.
    """
    # Replace the value and return the updated DataFrame
    return dataframe.replace(to_replace=target_value, value=replacement_value)


def replace_dataframe_contents(
    original_df: pd.DataFrame,
    new_df: pd.DataFrame,
    pad_missing: bool = False,
) -> pd.DataFrame:
    """
    Replace the contents of the original dataframe with the new dataframe.
    Optionally pad missing columns in the new dataframe and ensure all values are strings.

    :param original_df: The original dataframe with the expected headers.
    :param new_df: The new dataframe without headers.
    :param pad_missing: If True, pad missing columns in the new dataframe.
    :return: The updated dataframe with the same headers as the original.
    """
    # Ensure all values in new_df are converted to strings
    new_df = new_df.astype(str)

    # Clear the original dataframe
    original_headers = original_df.columns.tolist()

    if pad_missing:
        # Check for missing columns and pad them
        missing_columns = len(original_headers) - len(new_df.columns)
        if missing_columns > 0:
            logger.info(
                f"Padding {missing_columns} missing columns in the new dataframe."
            )
            for _ in range(missing_columns):
                new_df[len(new_df.columns)] = ""
        elif missing_columns < 0:
            raise ValueError(
                f"The new dataframe has more columns ({len(new_df.columns)}) "
                f"than the original dataframe ({len(original_headers)})."
            )

    # Validate column count
    if len(new_df.columns) != len(original_headers):
        raise ValueError(
            f"The new dataframe does not match the original dataframe's structure. "
            f"Expected {len(original_headers)} columns but got {len(new_df.columns)}."
        )

    # Assign headers from the original dataframe
    new_df.columns = original_headers

    # Return updated dataframe
    return new_df


def rename_data_column_content(
    section: Section,
    column_name: str,
    new_content: str,
    data_handler: DataHandler = DataHandler,
):
    data_handler = data_handler()
    data_handler.section = section
    data_handler.process(section)
    data_df = data_handler.content
    data_df[column_name] = new_content
    data_handler.content = data_df
    return data_handler.export()


def create_includes_section(
    include_path: str, include_handler: IncludesHandler = IncludesHandler
):
    include_handler = include_handler()
    # Create a dummy Section to initialize the handler
    dummy_section = Section(
        construct_name="include", handler_name="IncludesHandler", name=None
    )
    include_handler.section = dummy_section  # Assign the dummy Section

    # Set the content for the include line
    include_handler.content = f'#include "{include_path}"'

    # Export the section
    include_section = include_handler.export()
    return include_section


def get_gro_handler(
    gro_file: str,
    parser: GromacsParser = GromacsParser(),
    gro_handler: GroHandler = GroHandler(),
):
    sections = parser.parse(gro_file)
    gro_section = next(iter(sections.values()))
    gro_handler.process(gro_section)
    return gro_handler


def rename_residue_name_from_handler(gro_handler: GroHandler, new_residue_name: str):
    content_df = gro_handler.content
    content_df["Residue Name"] = new_residue_name

    # Update the handler content with the modified DataFrame
    gro_handler.content = content_df
    return gro_handler


def get_residue_number(gro_handler: GroHandler):
    residue_numbers = gro_handler.content["Residue Number"].iloc[-1]

    return residue_numbers


def validate_and_extract_residue_name(
    gro_handler: GroHandler, column_name="Residue Name"
):

    unique_values = gro_handler.content[column_name].unique()

    if len(unique_values) == 1:
        # Return the single unique value
        return unique_values[0]
    else:
        raise ValueError(
            f"Column '{column_name}' contains multiple residue names: {unique_values}"
        )


def rename_residue_name_from_gro(
    gro_file: str,
    new_residue_name: str,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
    gro_handler: GroHandler = GroHandler(),
):
    gro_handler = get_gro_handler(gro_file)
    gro_handler = rename_residue_name_from_handler(gro_handler, new_residue_name)

    output_file_path = prepare_output_file_path(
        gro_file, "gro", output_dir, output_name
    )
    output_file_path = export_gro_handler(gro_handler, output_file_path, parser)
    return output_file_path


def export_gro_handler(
    gro_handler: GroHandler,
    output_path: str,
    parser: GromacsParser = GromacsParser(),
):
    gro_section = gro_handler.export()
    sections = OrderedDict()
    sections["gro_file"] = gro_section
    output_file = parser.export(sections, output_path)
    return output_file


def add_full_rows_to_handler_deduplicate(
    handler: Union[GroHandler, DataHandler],
    dataframe_to_add: pd.DataFrame,
    add_to_top: bool = False,
    deduplicate_column: Optional[str] = None,
    keep: str = "first",
) -> Union[GroHandler, DataHandler]:
    """
    Add full rows to a handler's content, with optional deduplication based on a specified column.

    :param handler: The handler containing the original content.
    :param dataframe_to_add: The new dataframe to add to the handler's content.
    :param add_to_top: If True, adds the new rows to the top; otherwise, to the bottom.
    :param deduplicate_column: The column to use for deduplication. If None, no deduplication is performed.
    :param keep: Which duplicate to keep when deduplicating. "first" or "last".
    :return: The updated handler with the modified content.
    """
    # Combine the existing content with the new rows
    if add_to_top:
        combined_content = pd.concat(
            [dataframe_to_add, handler.content], ignore_index=True
        )
    else:
        combined_content = pd.concat(
            [handler.content, dataframe_to_add], ignore_index=True
        )

    # Deduplicate if a column is specified
    if deduplicate_column is not None:
        if deduplicate_column not in combined_content.columns:
            raise ValueError(
                f"Column '{deduplicate_column}' not found in the dataframe."
            )

        combined_content = combined_content.drop_duplicates(
            subset=deduplicate_column, keep=keep
        )

    # Update the handler content
    handler.content = combined_content
    return handler


def add_to_specific_handler_columns(
    handler: Union["GroHandler", "DataHandler"],
    column_values: Dict[str, Union[str, int, float]],
    add_to_top: bool = False,
) -> Union["GroHandler", "DataHandler"]:
    """
    Add a row with values for specific columns, leaving others empty.

    :param handler: The handler object containing a `content` DataFrame.
    :param column_values: Dictionary of column names and their values.
    :param add_to_top: Whether to add the row to the top or bottom of the DataFrame.
    :return: The updated handler.
    """
    # Create a row with specified column values, leaving others as NaN

    column_values_str = {key: str(value) for key, value in column_values.items()}

    # Create a row with specified column values, leaving others as NaN
    row_to_add = {
        col: column_values_str.get(col, None) for col in handler.content.columns
    }
    row_df = pd.DataFrame([row_to_add], columns=handler.content.columns)

    if add_to_top:
        new_content = pd.concat([row_df, handler.content], ignore_index=True)
    else:
        new_content = pd.concat([handler.content, row_df], ignore_index=True)
    handler.content = new_content
    return handler
