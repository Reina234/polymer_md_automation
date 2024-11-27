import os
import shutil
import logging
from typing import List

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
):
    """
    Copies specified files from source_dir to dest_dir. Optionally deletes the originals.

    Args:
        files_to_move (List[str]): List of specific file names to copy.
        source_dir (str): The directory to copy files from.
        dest_dir (str): The directory to copy files to.
        delete_original (bool): If True, deletes the files from the source_dir after copying.
    """
    # Ensure destination directory exists
    os.makedirs(dest_dir, exist_ok=True)

    for file_name in files_to_move:
        source_file = os.path.join(source_dir, file_name)
        dest_file = os.path.join(dest_dir, file_name)

        if os.path.exists(source_file):
            # Copy the file
            shutil.copy2(source_file, dest_file)
            logger.info(f"Copied {source_file} to {dest_file}.")

            # Delete the file if the flag is set
            if delete_original:
                os.remove(source_file)
                logger.info(f"Deleted original file: {source_file}")
        else:
            logger.warning(f"File not found: {source_file}")
