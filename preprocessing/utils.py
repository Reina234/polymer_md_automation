import os
import shutil
import logging
from typing import List
from pathlib import Path

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

logging.basicConfig(level=logging.ERROR)


def check_file_exists(input_file_path: str, get_contents: bool = False):
    if not os.path.exists(input_file_path):
        raise FileNotFoundError(f"Input file not found: {input_file_path}")

    with open(input_file_path, "r") as file:
        lines = file.readlines()

    if not lines:
        raise ValueError(f"PDB file is empty: {input_file_path}")

    if not get_contents:
        return
    else:
        return lines


def check_file_type(input_file_path: str, expected_file_type: str):
    """
    Validates that the file types of the input match
    the expected input_file_type

    Raises:
        ValueError: If validation fails.
    """

    input_extension = os.path.splitext(input_file_path)[1].lstrip(".")

    if input_extension != expected_file_type:
        message = (
            f"Validation failed: Expected input file of type "
            f"'{expected_file_type}', but got '{input_extension}'."
        )
        logging.error(message)
        raise ValueError(message)

    # Validation passed
    logging.info(f"Validation passed: Input file is of type '.{expected_file_type}' ")


def copy_files(
    files_to_move: List[str],
    source_dir: str,
    dest_dir: str,
    delete_original: bool = False,
) -> List[str]:
    """
    Copies specified files from source_dir to dest_dir. Optionally deletes the originals.

    Args:
        files_to_move (List[str]): List of specific file names to copy.
        source_dir (str): The directory to copy files from.
        dest_dir (str): The directory to copy files to.
        delete_original (bool): If True, deletes the files from the source_dir after copying.

    Returns:
        List[str]: List of file paths in the destination directory.
    """
    # Ensure destination directory exists
    os.makedirs(dest_dir, exist_ok=True)

    copied_files = []

    for file_name in files_to_move:
        source_file = os.path.join(source_dir, file_name)
        dest_file = os.path.join(dest_dir, file_name)

        if os.path.exists(source_file):
            # Copy the file
            shutil.copy2(source_file, dest_file)
            logger.info(f"Copied {source_file} to {dest_file}.")
            copied_files.append(dest_file)

            # Delete the file if the flag is set
            if delete_original:
                os.remove(source_file)
                logger.info(f"Deleted original file: {source_file}")
        else:
            logger.warning(f"File not found: {source_file}")

    return copied_files


import os
from pathlib import Path
from typing import List


def rename_basenames(file_paths: List[str], new_basename: str) -> List[str]:
    """
    Renames multiple files to have the same basename, preserving directories and extensions.

    Args:
        file_paths (List[str]): List of file paths to rename.
        new_basename (str): The new basename to apply to all files (without extensions).

    Returns:
        List[str]: List of renamed file paths.
    """
    renamed_files = []

    for file_path in file_paths:
        # Ensure the file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        # Get directory and extension
        directory = os.path.dirname(file_path)
        extension = os.path.splitext(file_path)[1]

        # Create the new file path
        new_file_path = os.path.join(directory, f"{new_basename}{extension}")

        # Skip renaming if already named correctly
        if os.path.abspath(file_path) == os.path.abspath(new_file_path):
            print(f"Skipping rename: '{file_path}' already has the desired basename.")
            renamed_files.append(new_file_path)
            continue

        # Rename the file
        Path(file_path).rename(new_file_path)
        renamed_files.append(new_file_path)

    return renamed_files


def OLD_rename_basenames(file_paths: List[str], new_basename: str) -> List[str]:
    """
    Renames multiple files to have the same basename with unique suffixes, preserving directories and extensions.

    Args:
        file_paths (List[str]): List of file paths to rename.
        new_basename (str): The new basename to apply to all files (without extensions).

    Returns:
        List[str]: List of renamed file paths.
    """
    renamed_files = []

    for i, file_path in enumerate(file_paths):
        # Ensure the file exists
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        # Get directory and extension
        directory = os.path.dirname(file_path)
        extension = os.path.splitext(file_path)[1]

        # Create the new file path with a unique suffix
        new_file_path = os.path.join(directory, f"{new_basename}_{i+1}{extension}")

        # Skip renaming if already named correctly
        if os.path.abspath(file_path) == os.path.abspath(new_file_path):
            print(f"Skipping rename: '{file_path}' already has the desired basename.")
            renamed_files.append(new_file_path)
            continue

        # Rename the file
        Path(file_path).rename(new_file_path)
        renamed_files.append(new_file_path)

    return renamed_files
