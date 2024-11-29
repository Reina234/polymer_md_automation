from config.paths import TOPOL_NAME, BASE_OUTPUT_DIR, GROMACS_OUTPUT_SUBDIR
import os
from pathlib import Path
from typing import Optional
from preprocessing.parsers.top_parser import TOPParser
from preprocessing.parsers.itp_parser import ITPParser
from preprocessing.utils import copy_files
from preprocessing.solvent_atomtypes_manager import AtomtypesManager
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def prepare_topol_file(
    input_top_path: str, run_name: str, output_base_dir: str = BASE_OUTPUT_DIR
):
    output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
    os.makedirs(output_dir, exist_ok=True)
    topol_file = move_and_rename_topol_file(input_top_path, output_dir, TOPOL_NAME)
    return topol_file


def move_and_rename_topol_file(
    input_top_path: str, output_dir: str, output_top_name: str
) -> str:
    """
    Copies and renames the topology file to the specified directory.

    Args:
        input_top_path (str): Path to the source topology file.
        output_dir (str): Path to the output directory.
        output_top_name (str): New name for the topology file.

    Returns:
        str: Full path to the renamed topology file.
    """
    # Ensure the source file exists
    if not os.path.exists(input_top_path):
        raise FileNotFoundError(
            f"The input topology file '{input_top_path}' does not exist."
        )

    top_name = os.path.basename(input_top_path)
    top_directory = os.path.dirname(input_top_path)

    # Copy the file to the output directory
    copied_files = copy_files([top_name], top_directory, output_dir)

    # Ensure at least one file is copied
    if not copied_files:
        raise ValueError(f"File '{top_name}' could not be copied to '{output_dir}'.")

    copied_file_path = Path(
        copied_files[0]
    )  # Extract the first file path from the list

    # Skip renaming if the source and destination paths are the same
    dest_file = os.path.join(output_dir, output_top_name)
    if os.path.abspath(copied_file_path) == os.path.abspath(dest_file):
        logger.info(
            f"Skipping rename: '{copied_file_path}' and '{dest_file}' are the same file."
        )
        return str(dest_file)

    # Rename the copied file
    renamed_file = copied_file_path.rename(dest_file)
    logger.info(f"Renamed '{copied_file_path}' to '{renamed_file}'.")
    return str(renamed_file)


def reformat_topol_file(
    input_top_path: str,
    solute_itp_path: str,
    solvent_itp_path: str,
    forcefield_path: str,
    output_dir: Optional[str] = None,
    posres: bool = False,
) -> str:
    """
    Reformat the topology file to ensure correct include order and handle position restraints.

    Args:
        input_top_path (str): Path to the input topology file.
        solute_itp_path (str): Path to the solute .itp file.
        solvent_itp_path (str): Path to the solvent .itp file.
        forcefield_path (str): Path to the forcefield .itp file.
        output_dir (Optional[str]): Directory to save the reformatted topology file.
        posres (bool): Whether to include position restraints.

    Returns:
        str: Path to the reformatted topology file.
    """
    # Ensure all paths are strings and resolve absolute paths
    if isinstance(solute_itp_path, list) or isinstance(solvent_itp_path, list):
        raise TypeError(
            "solute_itp_path and solvent_itp_path must be strings, not lists."
        )
    solute_itp_path = os.path.abspath(solute_itp_path)
    solvent_itp_path = os.path.abspath(solvent_itp_path)

    parser = TOPParser()

    # Ensure the input topology file is read
    content = parser.read_file(input_top_path)

    # Reformat includes and remove the defaults section
    content = parser.ensure_include_order(
        content,
        forcefield_include=forcefield_path,
        monomer_include=solute_itp_path,
        solvent_include=solvent_itp_path,
    )
    content = parser.handle_posres(content, posres)
    content = parser.remove_section(content, "defaults")

    # Save the updated topology file
    output_top_path = (
        os.path.join(output_dir, TOPOL_NAME) if output_dir else input_top_path
    )
    parser.save(output_top_path, content)
    return output_top_path


def add_atomtypes_to_topology(
    solvent_name: str, target_topology: str, include_header: bool = False
):
    """
    Add solvent atomtypes directly into the [ atomtypes ] section of the target topology file,
    while preserving any existing atomtypes.

    Args:
        solvent_name (str): Name of the solvent whose atomtypes to retrieve.
        target_topology (str): Path to the target topology file.
        include_header (bool): Whether to include the `[ atomtypes ]` header and comments.
    """
    manager = AtomtypesManager()
    new_atomtypes = manager.retrieve_atomtypes(solvent_name)

    # Optionally remove the `[ atomtypes ]` header and comments
    if not include_header:
        new_atomtypes = [
            line
            for line in new_atomtypes
            if not line.strip().startswith("[") and not line.strip().startswith(";")
        ]

    with open(target_topology, "r") as file:
        content = file.readlines()

    in_atomtypes_section = False
    updated_content = []
    existing_atomtypes = []
    atomtypes_added = False

    for line in content:
        stripped_line = line.strip()

        # Detect the start of the [ atomtypes ] section
        if stripped_line.startswith("[ atomtypes ]"):
            in_atomtypes_section = True
            updated_content.append(line)  # Keep the section header
            continue

        # Collect existing atomtypes within the section
        if in_atomtypes_section:
            if stripped_line.startswith("[") and not stripped_line.startswith(
                "[ atomtypes ]"
            ):
                in_atomtypes_section = False  # End of section
                updated_content.extend(existing_atomtypes)  # Add merged atomtypes
                updated_content.extend(new_atomtypes)  # Add new atomtypes
                atomtypes_added = True
                updated_content.append("\n")  # Add a blank line for separation
                updated_content.append(line)  # Keep the next section header
            else:
                existing_atomtypes.append(line)
                continue

        # Add lines to the updated content
        if not in_atomtypes_section:
            updated_content.append(line)

    # If no [ atomtypes ] section was found, create one at the beginning
    if not atomtypes_added:
        updated_content.insert(0, "\n")
        updated_content.insert(0, "\n".join(new_atomtypes) + "\n")
        updated_content.insert(0, "\n".join(existing_atomtypes) + "\n")
        updated_content.insert(0, "[ atomtypes ]\n")
        logger.info("[+] Created new [ atomtypes ] section in the topology file.")

    # Save the updated topology file
    with open(target_topology, "w") as file:
        file.writelines(updated_content)

    logger.info(
        f"[+] Added {solvent_name} atomtypes to the [ atomtypes ] section of {target_topology}."
    )
    return target_topology
