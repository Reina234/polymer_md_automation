import os
import shutil
import logging
from typing import List
from pathlib import Path
from functools import wraps
from typing import Optional
import logging

logger = logging.getLogger(__name__)


def check_file_exists(file_path: str, suppress_error: bool = False) -> Optional[str]:
    """
    Checks that a file exists in the correct location

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param suppress_error: If True, suppresses the FileNotFound error when file doesn't exist, and logs a warning instead, function will return None, defaults to False
    :type suppress_error: bool
    :raises FileNotFoundError: If the file is not found
    :return: The file path, if found, else None
    :rtype: Optional[str]
    """

    if not os.path.exists(file_path):
        message = f"File not found: {file_path}"
        if suppress_error:
            logger.warning(message)
            return None
        logger.error(message)
        raise FileNotFoundError(message)

    logger.info(f"File found: {file_path}")
    return file_path


def check_directory_exists(
    directory_path: str, make_dirs: bool = True
) -> Optional[str]:
    """
    Checks that a directory exists in the correct location

    :param directory_path: Relative or absolute path to the directory
    :type directory_path: str
    :param make_dirs: If True, missing directory is created, defaults to True
    :type make_dirs: bool, optional
    :raises FileNotFoundError: If the directory is not found
    :return: The directory path, if found, else None
    :rtype: Optional[str]
    """
    if not os.path.exists(directory_path):
        if make_dirs:
            logger.info(f"Creating directory: {directory_path}")

            os.makedirs(directory_path, exist_ok=True)

            logger.info(f"Directory created: {directory_path}")
            return directory_path
        else:
            logger.error(f"Directory not found: {directory_path}")

            raise FileNotFoundError(f"Directory not found: {directory_path}")
    else:
        logger.info(f"Directory found: {directory_path}")

        return directory_path


def check_file_contents(
    file_path: str, get_contents: bool = False
) -> Optional[List[str]]:
    """
    Checks that the contents of a file is not empty, and retrieves contents if :param get_contents: is True

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param get_contents: If True, returns the contents of the file as a list of strings, defaults to False
    :type get_contents: bool, optional
    :raises ValueError: If the file is empty
    :return: The contents of the file, if found and :param get_contents: is True, else None
    :rtype: Optional[List[str]]
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    if not lines:
        logger.error(f"File is empty: {file_path}")
        raise ValueError(f"File is empty: {file_path}")

    logger.info(f"File content checked: {file_path} is not empty")
    if not get_contents:
        return None
    else:
        return lines


def check_file_type(file_path: str, expected_file_type: str) -> None:
    """
    Validates that the file types of the input match the expected input_file_type

    :param file_path: Relative or absolute path to the file
    :type file_path: str
    :param expected_file_type: The expected file type
    :type expected_file_type: str
    :raises ValueError: If the file type does not match the expected file type
    """
    observed_file_type = os.path.splitext(file_path)[1].lstrip(".")

    if observed_file_type != expected_file_type:
        message = (
            f"Validation failed: Expected input file of type "
            f"'{expected_file_type}', but got '{observed_file_type}'."
        )
        logger.error(message)
        raise ValueError(message)

    logger.info(f"Validation passed: Input file is of type '.{expected_file_type}'")


def batch_copy_file(
    files_to_move: List[Optional[str]],
    dest_dir: str,
    delete_original: bool = False,
) -> List[Optional[str]]:
    """
    Copies specified files to a destination directory. Optionally deletes the originals.

    :param files_to_move: List of specific file names to copy, file name can be None, but will return None
    :type files_to_move: List[Optional[str]]
    :param dest_dir: The directory to copy files to
    :type dest_dir: str
    :param delete_original: If True, deletes the original files after copying, defaults to False
    :type delete_original: bool, optional
    :return: List of file paths in the destination directory
    :rtype: List[Optional[str]]
    """

    copied_files = []

    for file in files_to_move:
        copied_file = copy_file(
            file_path=file, dest_dir=dest_dir, delete_original=delete_original
        )
        copied_files.append(copied_file)

    return copied_files


def copy_file(
    file_path: Optional[str], dest_dir: str, delete_original: bool = False
) -> Optional[str]:
    """
    Copies a file to a destination directory. Optionally deletes the original.

    :param file_path: Relative or absolute path to the file to copy, can be None (will return None)
    :type file_path: Optional[str]
    :param dest_dir: The directory to copy the file to
    :type dest_dir: str
    :param delete_original: If True, deletes the original file after copying, defaults to False, defaults to False
    :type delete_original: bool, optional
    :return: The path to the copied file, if successful, else None
    :rtype: List[Optional[str,]]
    """

    if file_path is None:
        logger.info("File is None, skipping.")
        return None

    check_file_exists(file_path)
    check_directory_exists(dest_dir)

    dest_file = shutil.copy2(file_path, dest_dir)
    logger.info(f"Copied {file_path} to {dest_file}.")

    if delete_original:
        os.remove(file_path)
        logger.info(f"Deleted original file: {file_path}")

    return dest_file


def batch_rename_to_same(file_paths: List[str], new_name: str) -> List[Optional[str]]:
    """
    Renames multiple files to have the same name, preserving directories and extensions.

    :param file_paths: List of specific file names to rename, file name can be None, but will return None
    :type file_paths: List[str]
    :param new_name: The new name for ALL the files, if you want to specify names, please use batch_rename_to_list()
    :type new_name: str
    :return: List of renamed file paths
    :rtype: List[Optional[str]]
    """
    renamed_files = []

    for file in file_paths:
        renamed_file = rename_file(file_path=file, new_name=new_name)
        renamed_files.append(renamed_file)

    return renamed_files


def batch_rename_to_list(
    file_paths: List[str], new_names: List[str]
) -> List[Optional[str]]:
    """
    Renames multiple files to have the names specified in new_names, preserving directories and extensions.

    :param file_paths: List of specific file names to rename, file name can be None, but will return None
    :type file_paths: List[str]
    :param new_names: List of new names for the files, if you want to rename it all to the same name, please use batch_rename_to_same()
    :type new_names: List[str]
    :return: List of renamed file paths
    :rtype: List[Optional[str]]
    """

    renamed_files = []

    for new_name, file in zip(new_names, file_paths):
        renamed_file = rename_file(file_path=file, new_name=new_name)
        renamed_files.append(renamed_file)

    return renamed_files


def rename_file(file_path: Optional[str], new_name: str) -> Optional[str]:
    """
    Renames a file to a new name.

    :param file_path: Relative or absolute path to the file to rename, can be None (will return None)
    :type file_path: Optional[str]
    :param new_name: The new name for the file
    :type new_name: str
    :return: The new file path
    :rtype: str
    """
    if file_path is None:
        logger.info("File is None, skipping.")
        return None

    file_path = check_file_exists(file_path)
    path = Path(file_path)

    new_name_with_file_extension = new_name + path.suffix
    new_file_path = path.with_name(new_name_with_file_extension)

    if new_file_path.exists():
        raise FileExistsError(f"A file with the name '{new_file_path}' already exists.")

    os.rename(file_path, new_file_path)

    logger.info(f"Renamed {file_path} to {new_file_path}")
    return str(new_file_path)
