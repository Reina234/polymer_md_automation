from config.paths import TOPOL_NAME, BASE_OUTPUT_DIR, GROMACS_OUTPUT_SUBDIR
import os
from pathlib import Path
from typing import Optional
from preprocessing.parsers.top_parser import TOPParser
from preprocessing.parsers.itp_parser import ITPParser
from preprocessing.utils import copy_files
from preprocessing.solvent_atomtypes_manager import SolventAtomtypesManager


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
    copied_file_path = Path(copy_files([top_name], top_directory, output_dir))

    # Skip renaming if the source and destination paths are the same
    dest_file = os.path.join(output_dir, output_top_name)
    if os.path.abspath(copied_file_path) == os.path.abspath(dest_file):
        print(
            f"Skipping rename: '{copied_file_path}' and '{dest_file}' are the same file."
        )
        return str(dest_file)

    # Rename the copied file
    renamed_file = copied_file_path.rename(dest_file)
    return str(renamed_file)


def reformat_topol_file(
    input_top_path: str,
    solute_itp_path: str,
    solvent_itp_path: str,
    forcefield_path: str,
    output_dir: Optional[str] = None,
    posres: bool = False,
) -> str:
    parser = TOPParser()
    solute_itp_path = os.path.abspath(solute_itp_path)
    solvent_itp_path = os.path.abspath(solvent_itp_path)
    forcefield_include = forcefield_path
    solvent_include = solvent_itp_path
    solute_include = solute_itp_path
    content = parser.read_file(input_top_path)
    content = parser.ensure_include_order(
        content, forcefield_include, solute_include, solvent_include
    )
    parser.handle_posres(content, posres)
    content = parser.remove_section(content, "defaults")

    output_top_path = (
        os.path.join(output_dir, TOPOL_NAME) if output_dir else input_top_path
    )

    parser.save(output_top_path, content)
    return output_top_path


def add_atomtypes_to_topology(solvent_name: str, target_topology: str):
    manager = SolventAtomtypesManager()
    atomtypes = manager.retrieve_atomtypes(solvent_name)

    # Read target topology content
    with open(target_topology, "r") as file:
        target_content = file.readlines()

    # Append atomtypes to the target topology
    target_content.append("\n")  # Ensure spacing
    target_content.extend(atomtypes)

    # Save the updated topology
    with open(target_topology, "w") as file:
        file.writelines(target_content)
    logger.info(f"[+] Added {solvent_name} atomtypes to {target_topology}.")
