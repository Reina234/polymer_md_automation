from abc import ABC, abstractmethod
from typing import Optional, Tuple
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
import logging
from config.constants import LengthUnits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class BaseGromacsCommand(ABC):
    OUTPUT_NAME: str  # Subclasses must define this as a class variable
    UNITS = LengthUnits.NANOMETER

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def __init_subclass__(cls):
        super().__init_subclass__()
        if not hasattr(cls, "OUTPUT_NAME"):
            raise TypeError(f"{cls.__name__} must define a class-level 'OUTPUT_NAME'")

    @abstractmethod
    def run(self, **kwargs) -> str:
        pass

    def _execute(self, command) -> str:
        try:
            logger.info(f"Starting GROMACS command: {self.__class__.__name__}")
            logger.debug(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True)
            logger.info(f"Command completed successfully: {self.__class__.__name__}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running GROMACS command: {e}")
            raise

    @abstractmethod
    def metadata(self, **kwargs) -> dict:
        """
        Generates metadata for the GROMACS command execution.

        Returns:
            dict: Metadata details.
        """
        pass

    def _update_metadata(self, **kwargs):
        """
        Updates metadata using the provided dictionary.

        Args:
            metadata (dict): Metadata details to record.
        """
        metadata = self.metadata(**kwargs)
        self.metadata_tracker.add_step(
            step_name=self.__class__.__name__, details=metadata
        )
