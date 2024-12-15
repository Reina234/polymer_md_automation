from abc import ABC, abstractmethod
from typing import Optional, Tuple
import subprocess
from A_modules.shared.metadata_tracker import MetadataTracker
import logging
from config.constants import LengthUnits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class CommandLineOperation(ABC):
    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    @abstractmethod
    def run(self, **kwargs) -> str:
        pass

    def _execute(self, command) -> str:
        try:
            logger.info(f"Starting command: {self.__class__.__name__}")
            logger.debug(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True)
            logger.info(f"Command completed successfully: {self.__class__.__name__}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running  command: {e}")
            raise

    @abstractmethod
    def metadata(self, **kwargs) -> dict:
        """
        Generates metadata for the  command execution.

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
