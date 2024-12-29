from A_modules.shared.utils.dataframe_utils import (
    dataframe_not_empty_check,
)
from A_modules.shared.utils.file_utils import (
    file_exists_check_wrapper,
    file_type_check_wrapper,
    check_directory_exists,
    prepare_output_file_path,
    add_suffix_to_filename,
    copy_file,
)
from A_modules.atomistic.gromacs.parser.handlers.includes_handler import IncludesHandler
from collections import OrderedDict
from A_modules.atomistic.gromacs.commands.insert_molecules import InsertMolecules
from A_config.paths import TEMP_DIR
from A_config.constants import MassUnits2, LengthUnits2
from A_modules.shared.utils.calculation_utils import calculate_num_particles
import pandas as pd
from A_modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from A_modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler
from typing import List, Optional, Dict, Union
from data_models.solvent import Solvent
from A_modules.atomistic.gromacs.commands.solvate import Solvate
from A_modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
import logging
import os
from A_modules.atomistic.gromacs.commands.editconf import Editconf
from data_models.output_types import GromacsPaths
from A_modules.atomistic.gromacs.parser.data_models.section import Section

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
    residue_number_col="Residue Number",
) -> pd.DataFrame:
    """
    Calculates the number of molecules for each unique residue name by counting unique residue numbers.

    :param dataframe: Input dataframe containing residue information.
    :param residue_number_col: Column name for residue numbers.
    :param residue_name_col: Column name for residue names.
    :return: A dataframe with residue names and their corresponding molecule counts.
    """
    dataframe = gro_handler.content

    # Group by residue name and count unique residue numbers
    molecule_counts = (
        dataframe.groupby(residue_name_col)[residue_number_col]
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


@dataframe_not_empty_check(dataframe_arg_index=0)
def calculate_minimum_box_size_from_df(
    atom_df: pd.DataFrame,
    padding: float = 0.1,
    expected_columns: List[str] = ["X", "Y", "Z"],
) -> List[float]:
    if atom_df.empty:
        raise ValueError("Atom data is empty. Cannot calculate bounding box size.")

    if not all(col in atom_df.columns for col in expected_columns):
        raise ValueError(
            f"DataFrame must contain columns: {', '.join(expected_columns)}"
        )

    # Calculate min and max for each dimension
    min_coords = atom_df[expected_columns].min().values
    max_coords = atom_df[expected_columns].max().values

    # Add padding and calculate box size
    box_size = (max_coords - min_coords) + 2 * padding
    return box_size.tolist()


@file_exists_check_wrapper(file_arg_index=0)
def count_particles(gro_file: str) -> int:
    parser = GromacsParser()
    sections = parser.parse(gro_file)
    gro_section = next(iter(sections.values()))
    gro_handler = GroHandler()
    gro_handler.process(gro_section)
    residue_numbers = gro_handler.content["Residue Number"].unique()
    return len(residue_numbers)


def num_molecules_needed(gro_file: str, desired_number: int) -> int:
    current_particles = count_particles(gro_file)
    return desired_number - current_particles


# @file_type_check_wrapper(file_arg_index=0, expected_file_type="gro")
def calculate_minimum_box_size_from_gro(
    gro_file, padding: float = 0.1, parser=GromacsParser(), gro_handler=GroHandler()
) -> List[float]:
    sections = parser.parse(gro_file)
    gro_section = next(iter(sections.values()))
    gro_handler.process(gro_section)
    box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding)
    return box_size


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


# NOTE: honestly, could turn the grohandler -> dataframe, dataframe -> grohandler part into a wrapper or something
def add_full_rows_to_handler(
    handler: Union[GroHandler, DataHandler],
    dataframe_to_add: pd.DataFrame,
    add_to_top: bool = False,
) -> Union[GroHandler, DataHandler]:

    if add_to_top is True:
        new_content = pd.concat([dataframe_to_add, handler.content], ignore_index=True)
    else:
        new_content = pd.concat([handler.content, dataframe_to_add], ignore_index=True)
    handler.content = new_content
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


def calculate_minimum_box_size_from_handler(
    gro_handler: GroHandler, padding: float = 0.1
):
    box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding)
    return box_size


def check_or_create_box_dim(gro_handler: GroHandler):
    if not gro_handler.box_dimensions:
        logging.warning(
            "Box dimensions not found in GRO file. Calculating minimum box size."
        )
        box_size = calculate_minimum_box_size_from_handler(gro_handler)
    else:
        box_size = gro_handler.box_dimensions
    return box_size


def validate_box_dimensions(
    box_dimensions: Optional[List[float]], suppress_error: bool = True
) -> Optional[List[float]]:
    """
    Validates that box_dimensions is a list of three floats.
    Handles NoneType gracefully.

    :param box_dimensions: List of box dimensions (or None).
    :param suppress_error: If True, logs a warning instead of raising an error.
    :return: The box dimensions if valid, otherwise None.
    """
    # Handle NoneType input
    if box_dimensions is None:
        logger.warning("Box dimensions are None. Returning None.")
        return None

    # Validate box dimensions are a list of three floats
    if not isinstance(box_dimensions, list) or len(box_dimensions) != 3:
        if not suppress_error:
            raise ValueError("Box dimensions must be a list of three floats.")
        else:
            logger.warning("Box dimensions are not a list of three floats.")
            return None

    if not all(isinstance(dim, float) for dim in box_dimensions):
        if not suppress_error:
            raise ValueError("All box dimensions must be of type float.")
        else:
            logger.warning("Box dimensions contain non-float values.")
            return None

    logger.info("Box dimensions are valid.")
    return box_dimensions


def validate_solute_gro_with_editconf(
    gro_file: str,
    output_dir: Optional[str] = None,
    validated_file_name: Optional[str] = None,
    suppress_error: bool = True,
    editconf: Editconf = Editconf(),
) -> str:
    if not output_dir:
        output_dir = os.path.dirname(gro_file)
        print("not new")
        print(type(output_dir))
    else:
        check_directory_exists(output_dir, make_dirs=True)
        print("new output dir")
        print(output_dir)

    gro_handler = get_gro_handler(gro_file)
    box_size = gro_handler.box_dimensions
    if not box_size:

        box_size = validate_box_dimensions(
            gro_handler.box_dimensions, suppress_error=suppress_error
        )
    if not validated_file_name:
        validated_file_name = add_suffix_to_filename(gro_file, "_validated")

    if box_size:
        return gro_file
    else:
        box_size = calculate_minimum_box_size_from_handler(gro_handler)
        print("TEST AGAIN !!!")
        print(gro_file)
        print(output_dir)
        validated_solute_path = editconf.run(
            gro_file,
            output_dir,
            box_size_nm=box_size,
            output_name=validated_file_name,
        )
        return validated_solute_path


@file_type_check_wrapper(file_arg_index=0, expected_file_type="gro")
@file_type_check_wrapper(file_arg_index=1, expected_file_type="top")
def prepare_solvated_solute_box(
    solvent_gro_file: str,
    topol_file: str,
    final_box_size: List[float],
    output_dir: str,
    editconf: Editconf = Editconf(),
    solvate: Solvate = Solvate(),
    temp_dir: str = TEMP_DIR,
    initial_box_factor: float = 0.9,
) -> str:

    initial_box_size_nm = [dim * initial_box_factor for dim in final_box_size]
    editconf_box_gro = editconf.run(
        solvent_gro_file,
        temp_dir,
        box_size_nm=initial_box_size_nm,
    )

    print(output_dir)
    print("finished editconf")
    solvated_box = solvate.run(
        editconf_box_gro,
        solvent_gro_file,
        topol_file,
        output_dir,
    )
    print("finished solvate")
    return solvated_box


def count_particles(gro_handler: GroHandler) -> int:
    residue_numbers = gro_handler.content["Residue Number"].unique()
    return len(residue_numbers)


def get_remaining_particles(
    gro_handler: GroHandler, target_num_particles: int, suppress_error: bool
) -> int:
    current_particles = count_particles(gro_handler)
    remaining_particles = int(target_num_particles - current_particles)
    if remaining_particles < 0:
        if not suppress_error:
            raise ValueError("Number of remaining particles is negative.")
        else:
            logging.warning("Number of remaining particles is negative.")
            return None
    return remaining_particles


def is_target_achieved(
    remaining_particles: int, target_num_particles: int, tolerance: float
):
    """
    Checks if the target number of particles has been achieved.
    """
    return remaining_particles <= target_num_particles * tolerance


@file_type_check_wrapper(file_arg_index=0, expected_file_type="gro")
def validated_solvated_box(solvent_box_gro: str, target_num_particles: int):
    gro_handler = get_gro_handler(solvent_box_gro)
    remaining_particles = get_remaining_particles(
        gro_handler, target_num_particles, suppress_error=False
    )
    if remaining_particles is None:
        return False


def scale_box_to_final_size(
    gro_file: str,
    final_box_size: List[float],
    output_dir: str,
    editconf: Editconf = Editconf(),
    output_name: Optional[str] = None,
) -> str:
    """
    Scales the box to the final size using editconf.

    :param gro_file: Path to the .gro file to scale.
    :param final_box_size: Desired final box size in nanometers.
    :param output_dir: Directory to save the output file.
    :param editconf: Editconf object for running the command.
    :param output_name: Optional name for the scaled .gro file.
    :return: Path to the scaled .gro file.
    """
    if not output_name:
        output_name = add_suffix_to_filename(gro_file, "scaled", return_full_path=False)

    logger.info(f"Scaling box to final size {final_box_size} nm.")
    scaled_box_path = editconf.run(
        gro_file,
        output_dir,
        box_size_nm=final_box_size,
        output_name=output_name,
    )
    return scaled_box_path


def create_initial_valid_box(
    solvent_gro_file: str,
    topol_file: str,
    final_box_size: List[float],
    target_num_particles: int,
    output_dir: str,
    initial_box_factor: float = 0.9,
    safety_margin: float = 0.95,
    max_attempts: int = 5,
    editconf: Editconf = Editconf(),
    solvate: Solvate = Solvate(),
) -> str:
    """
    Creates an initial solvated box that does not exceed the target particle count.

    :return: Path to the valid initial solvated box.
    """
    box_factor = initial_box_factor
    attempt = 0

    while attempt < max_attempts:
        attempt += 1
        logger.info(
            f"Attempt {attempt}: Creating solvated box with factor {box_factor:.3f}."
        )

        # Prepare the solvated box
        solvated_box_path = prepare_solvated_solute_box(
            solvent_gro_file,
            topol_file,
            final_box_size=final_box_size,
            output_dir=output_dir,
            initial_box_factor=box_factor,
            editconf=editconf,
            solvate=solvate,
        )

        # Validate particle count
        current_particles = count_particles(get_gro_handler(solvated_box_path))
        if current_particles <= target_num_particles:
            logger.info("Valid box created with acceptable particle count.")

            # Scale the box back to its full size
            scaled_box_path = scale_box_to_final_size(
                gro_file=solvated_box_path,
                final_box_size=final_box_size,
                output_dir=output_dir,
                editconf=editconf,
            )
            return scaled_box_path

        # Reduce box size factor and retry
        logger.warning(
            f"Particle count exceeded target ({current_particles} > {target_num_particles}). "
            "Reducing box size factor."
        )
        box_factor *= safety_margin

    raise RuntimeError(
        f"Failed to create a valid solvated box within {max_attempts} attempts."
    )


def add_missing_molecules(
    box_gro_file: str,
    solvent_gro_file: str,
    target_num_particles: int,
    tolerance: float,  # 1% tolerance by default
    max_attempts: int = 10,
    insert_molecules: InsertMolecules = InsertMolecules(),
) -> str:
    """
    Adds remaining molecules iteratively until the target particle count is reached within tolerance.

    :param box_gro_file: Path to the box .gro file.
    :param solvent_gro_file: Path to the solvent .gro file.
    :param target_num_particles: Desired number of particles.
    :param tolerance: Tolerance for the number of particles (default: 1%).
    :param max_attempts: Maximum number of attempts to add molecules.
    :return: Path to the updated box .gro file.
    """
    attempt = 0

    while attempt < max_attempts:
        attempt += 1
        logger.info(f"Attempt {attempt}: Checking remaining particles...")

        gro_handler = get_gro_handler(box_gro_file)
        remaining_particles = get_remaining_particles(
            gro_handler, target_num_particles, suppress_error=False
        )

        if is_target_achieved(remaining_particles, target_num_particles, tolerance):
            logger.info("Target number of particles achieved within tolerance.")
            return box_gro_file

        if remaining_particles > 0:
            logger.info(f"Adding {remaining_particles} molecules to the box.")
            box_gro_file = insert_molecules.run(
                box_gro_file,
                solvent_gro_file,
                num_molecules=remaining_particles,
            )
        else:
            logger.warning("No remaining particles to add. Stopping early.")
            break

    raise RuntimeError(
        f"Failed to achieve target particle count within {max_attempts} attempts."
    )


def refine_solvated_box(
    solvent_gro_file: str,
    topol_file: str,
    final_box_size: List[float],
    target_num_particles: int,
    output_dir: str,
    initial_box_factor: float,
    tolerance: float,
    safety_margin: float,
    max_attempts: int,
) -> str:
    """
    Refines a solvated box to meet the target particle count.

    :param solvent_gro_file: Path to the solvent .gro file.
    :param topol_file: Path to the topology file.
    :param final_box_size: Desired final box size in nanometers.
    :param target_num_particles: Target number of particles.
    :param output_dir: Directory to save the output.
    :param initial_box_factor: Initial scaling factor for the box.
    :param safety_margin: Safety margin for box size reduction.
    :param max_attempts: Maximum attempts to refine the box.
    :return: Path to the refined solvated box.
    """
    # Step 1: Create an initial valid box
    logger.info("Starting box refinement process...")
    box_gro_path = create_initial_valid_box(
        solvent_gro_file=solvent_gro_file,
        topol_file=topol_file,
        final_box_size=final_box_size,
        target_num_particles=target_num_particles,
        output_dir=output_dir,
        initial_box_factor=initial_box_factor,
        safety_margin=safety_margin,
        max_attempts=max_attempts,
    )

    # Step 2: Add remaining molecules if needed
    box_gro_path = add_missing_molecules(
        box_gro_file=box_gro_path,
        solvent_gro_file=solvent_gro_file,
        target_num_particles=target_num_particles,
        tolerance=tolerance,
    )

    # Step 3: Return the final refined box
    logger.info(f"Box refinement complete. Final box saved at: {box_gro_path}")
    return box_gro_path


def create_solvated_box(
    solvent_gro_file: str,
    topol_file: str,
    final_box_size: List[float],
    solvent: Solvent,
    output_dir: str,
    initial_box_factor: float = 0.9,
    tolerance: float = 0.01,
    safety_margin: float = 0.95,
    max_attempts: int = 5,
) -> str:
    validated_solvent_gro_file = validate_solute_gro_with_editconf(
        gro_file=solvent_gro_file, output_dir=output_dir
    )
    target_num_particles = calculate_num_particles(
        box_dimensions=final_box_size,
        molecular_weight=solvent.molecular_weight,
        density_SI=solvent.density,
        box_units=LengthUnits2.NANOMETER,
        mass_units=MassUnits2.GRAM,
    )
    solvated_box = refine_solvated_box(
        solvent_gro_file=validated_solvent_gro_file,
        topol_file=topol_file,
        final_box_size=final_box_size,
        target_num_particles=target_num_particles,
        output_dir=output_dir,
        initial_box_factor=initial_box_factor,
        tolerance=tolerance,
        safety_margin=safety_margin,
        max_attempts=max_attempts,
    )

    print("NUM_MOLECULES")
    return solvated_box


def prepare_solvent_box_name(solvent: Solvent, extension: str):
    return f"{solvent.name.lower()}_solvent_box.{extension}"


from A_modules.shared.file_conversion.converters.base_converter import BaseConverter
from A_modules.shared.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)
from A_modules.shared.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
from A_modules.shared.packmol.solvent_box import PackmolSolventBox


def create_solvent_box_gro(
    input_gro_file,
    output_dir: str,
    box_size_nm: List[float],
    solvent: Solvent,
    output_name: Optional[str] = None,
    temp_dir: str = TEMP_DIR,
    gro_to_pdb_converter: BaseConverter = EditconfGROtoPDBConverter(),
    pdb_to_gro_converter: BaseConverter = EditconfPDBtoGROConverter(),
    packmol_operation: PackmolSolventBox = PackmolSolventBox(),
) -> str:
    output_pdb = gro_to_pdb_converter.run(input_gro_file, temp_dir)
    if not output_name:
        output_name = prepare_solvent_box_name(solvent, "gro")

    packmol_output = packmol_operation.run(
        output_pdb,
        output_dir=temp_dir,
        solvent=solvent,
        box_size_nm=box_size_nm,
    )
    print(packmol_output)
    output_gro = pdb_to_gro_converter.run(
        packmol_output, output_dir, box_size_nm=box_size_nm, output_name=output_name
    )
    return output_gro
