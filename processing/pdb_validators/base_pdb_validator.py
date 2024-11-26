import logging
from abc import ABC, abstractmethod
from typing import Optional
from processing.metadata_tracker import MetadataTracker

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class BasePDBValidator(ABC):
    """
    Abstract base class for PDB file validators.
    """

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def add_metadata_tracker(self, metadata_tracker: MetadataTracker) -> None:
        """
        Attach a metadata tracker to the validator.
        """
        self.metadata_tracker = metadata_tracker

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

    def _skeletal_check(self, file_path: str) -> bool:
        """
        Perform the minimum validation required for a standard PDB file.
        Args:
            file_path (str): Path to the PDB file.

        Returns:
            bool: True if skeletal validation passes, False otherwise.
        """
        try:
            with open(file_path, "r") as file:
                lines = file.readlines()

            # Check for ATOM/HETATM records
            if not any(line.startswith(("ATOM", "HETATM")) for line in lines):
                logger.error(f"[!] PDB file {file_path} lacks ATOM or HETATM records.")
                return False

            # Check for END or ENDMDL
            if lines[-1].strip() not in {"END", "ENDMDL"}:
                logger.error(f"[!] PDB file {file_path} must end with 'END' or 'ENDMDL'.")
                return False

            logger.info(f"[+] Skeletal validation passed for {file_path}.")
            return True

        except Exception as e:
            logger.exception(f"[!] Skeletal validation failed for {file_path}: {e}")
            return False

    @property
    def supports_fixing(self) -> bool:
        """
        Indicates whether this validator supports fixing.
        """
        return False

    def _make_valid(self, file_path: str) -> None:
        """
        Attempt to fix the PDB file to pass validation.
        Raises:
            NotImplementedError: If fixing is not supported.
        """
        if not self.supports_fixing:
            raise NotImplementedError(f"[!] Fixing is not supported for {self.__class__.__name__}.")
        raise NotImplementedError(f"[!] `make_valid` must be implemented for {self.__class__.__name__}.")
