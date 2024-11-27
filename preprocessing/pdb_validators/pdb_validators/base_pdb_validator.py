import os
import logging
from typing import List, Optional
from preprocessing.file_commenter import FileCommenter
from preprocessing.utils import check_file_exists, check_file_type
from abc import ABC, abstractmethod

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class BasePDBValidator(ABC):
    """
    Performs a skeletal validation of a PDB file.
    Ensures the file has a valid format minimally compatible with tools like Packmol or PyMOL.
    """

    def __init__(self, file_commenter: Optional[FileCommenter] = None):
        self.file_commenter = file_commenter or FileCommenter()

    @abstractmethod
    def validate(self, file_path: str) -> bool:
        """
        Validate the PDB file.

        Args:
            file_path (str): Path to the PDB file.

        Returns:
            bool: True if the file passes validation, False otherwise.
        """
        pass

    def _skeletal_check(self, file_path: str) -> None:
        """
        Perform skeletal validation of a PDB file.

        Args:
            file_path (str): Path to the PDB file.

        Raises:
            ValueError: If the file fails skeletal validation.
        """
        try:
            lines = self._read_file(file_path)

            if not self._has_valid_atoms(lines):
                logger.error(f"[!] Validation failed: No valid ATOM or HETATM records found in {file_path}.")
                raise ValueError(f"File {file_path} failed skeletal validation.")

            logger.info(f"[+] Validation passed for {file_path}.")
            self.file_commenter.add_comment(file_path, "Validation passed: Basic skeletal check.")
        except Exception as e:
            logger.error(f"[!] Skeletal validation failed for {file_path}: {e}")
            raise  # Re-raise the exception to propagate it
        
        
    def _read_file(self, input_file_path: str) -> List[str]:
        """
        Read the PDB file into memory.

        Args:
            file_path (str): Path to the PDB file.

        Returns:
            List[str]: List of lines in the file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file cannot be read.
        """
        check_file_exists(input_file_path)
        lines = check_file_exists(input_file_path, "pdb")
        return lines

    def _has_valid_atoms(self, lines: List[str]) -> bool:
        """
        Check for valid ATOM or HETATM records in the file.

        Args:
            lines (List[str]): Lines from the PDB file.

        Returns:
            bool: True if at least one valid ATOM or HETATM record exists, False otherwise.
        """
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    float(line[30:38].strip())  # x-coordinate
                    float(line[38:46].strip())  # y-coordinate
                    float(line[46:54].strip())  # z-coordinate
                    return True
                except ValueError:
                    logger.warning(f"[!] Invalid atom coordinates in line: {line.strip()}")
        return False
