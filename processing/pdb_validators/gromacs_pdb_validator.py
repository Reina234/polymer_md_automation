from typing import Optional
from processing.pdb_validators.base_pdb_validator import BasePDBValidator
from processing.pdb_validators.pdb_utils import add_box_information, calculate_box_size
from data_models.solvent import Solvent
from processing.metadata_tracker import MetadataTracker
import logging

logger = logging.getLogger(__name__)

class GROMACSPDBValidator(BasePDBValidator):
    """
    Validator for GROMACS-compatible PDB files, currently used for solvation box validation.
    """

    def __init__(self, solvent: Solvent, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)
        self.solvent = solvent

    def validate(self, file_path: str) -> bool:
        """
        Validate the PDB file for GROMACS compatibility.
        Args:
            file_path (str): Path to the PDB file.

        Returns:
            bool: True if the file is valid, False otherwise.
        """
        if not self._skeletal_check(file_path):
            logger.error(f"[!] Skeletal validation failed for {file_path}.")
            return False

        if not self._has_box_information(file_path):
            logger.warning(f"[!] Box size information (CRYST1) missing in {file_path}.")
            if self.supports_fixing:
                try:
                    logger.info("[+] Attempting to fix the file.")
                    self._make_valid(file_path)
                    logger.info(f"[+] File fixed successfully: {file_path}.")
                    return True
                except Exception as e:
                    logger.error(f"[!] Failed to fix the file: {e}")
                    return False
            else:
                logger.error("[!] Fixing not supported for this validator.")
                return False

        logger.info(f"[+] GROMACS validation passed for {file_path}.")
        return True

    @property
    def supports_fixing(self) -> bool:
        """
        Indicates that this validator supports fixing.
        """
        return True

    def _make_valid(self, file_path: str) -> None:
        """
        Fix the PDB file for GROMACS compatibility using solvent properties.
        Args:
            file_path (str): Path to the PDB file.
        """
        if not self._has_box_information(file_path):
            box_size = calculate_box_size(self.solvent.molecular_weight, self.solvent.density)
            add_box_information(file_path, box_size)

        # Revalidate after fixing
        if not self.validate(file_path):
            raise RuntimeError(f"[!] PDB file {file_path} could not be fixed.")

    def _has_box_information(self, file_path: str) -> bool:
        """
        Check if the PDB file includes box information (CRYST1 line).
        Args:
            file_path (str): Path to the PDB file.

        Returns:
            bool: True if the CRYST1 line is present, False otherwise.
        """
        with open(file_path, "r") as file:
            return any(line.startswith("CRYST1") for line in file)
