import os
import logging
from typing import List, Optional
from preprocessing.file_commenter import FileCommenter
from preprocessing.utils import check_file_exists, check_file_type
from abc import ABC, abstractmethod
from data_models.solvent import Solvent
from preprocessing.pdb_validators.pdb_parser import PDBParser

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class BasePDBValidator(ABC):
    """
    Performs a skeletal validation of a PDB file.
    Ensures the file has a valid format minimally compatible with tools like Packmol or PyMOL.
    """

    SKELETAL_VALIDATION_PASSED_COMMENT = "Validation passed: Basic skeletal check."

    def __init__(self):
        self.pdb_parser = PDBParser()

    @abstractmethod
    def validate(
        self,
        input_file_path: str,
        solvent: Solvent,
        output_file_path: Optional[str] = None,
    ) -> bool:
        """
        Validate the PDB file.

        Args:
            file_path (str): Path to the PDB file.

        Returns:
            bool: True if the file passes validation, False otherwise.
        """
        pass

    def _skeletal_check(
        self, input_file_path: str, output_file_path: Optional[str] = None
    ) -> None:
        """
        Perform skeletal validation of a PDB file.

        Args:
            file_path (str): Path to the PDB file.

        Raises:
            ValueError: If the file fails skeletal validation.
        """
        if not output_file_path:
            output_file_path = input_file_path
        try:
            content = self.pdb_parser.read_file(input_file_path)

            if not self._has_valid_atoms(content):
                logger.error(
                    f"[!] Validation failed: No valid ATOM or HETATM records found in {input_file_path}."
                )
                raise ValueError(f"File {input_file_path} failed skeletal validation.")

            logger.info(f"[+] Validation passed for {input_file_path}.")
            content = self.pdb_parser.add_comment(
                content, self.SKELETAL_VALIDATION_PASSED_COMMENT
            )
            self.pdb_parser.save(output_file_path, content)

        except Exception as e:
            logger.error(f"[!] Skeletal validation failed for {input_file_path}: {e}")
            raise  # Re-raise the exception to propagate it

    def _has_valid_atoms(self, content: List[str]) -> bool:
        """
        Check for valid ATOM or HETATM records in the file.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            bool: True if at least one valid ATOM or HETATM record exists, False otherwise.
        """
        for line in content:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    float(line[30:38].strip())  # x-coordinate
                    float(line[38:46].strip())  # y-coordinate
                    float(line[46:54].strip())  # z-coordinate
                    return True
                except ValueError:
                    logger.warning(
                        f"[!] Invalid atom coordinates in line: {line.strip()}"
                    )
        return False
