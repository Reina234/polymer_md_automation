# File: utils/file_operations.py

import os
import shutil
import logging

logger = logging.getLogger(__name__)

def move_file(file_path: str, target_directory: str) -> str:
    """
    Move a file to a specified directory.

    Args:
        file_path (str): Path to the file to move.
        target_directory (str): Directory to move the file to.

    Returns:
        str: Path to the moved file.
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    os.makedirs(target_directory, exist_ok=True)
    target_path = os.path.join(target_directory, os.path.basename(file_path))
    shutil.move(file_path, target_path)
    logger.info(f"[+] File moved to {target_path}.")
    return target_path
