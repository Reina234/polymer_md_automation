from abc import ABC, abstractmethod
from preprocessing.metadata_tracker import MetadataTracker
from typing import Optional
import subprocess
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class GromacsCommand(ABC):
    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def __init_subclass__(cls):
        super().__init_subclass__()
        if not hasattr(cls, "OUTPUT_GRO_NAME"):
            raise TypeError(
                f"{cls.__name__} must define a class-level 'OUTPUT_GRO_NAME'"
            )

    def run(self, *args, **kwargs):
        try:
            logger.info(f"Starting GROMACS command: {self.__class__.__name__}")
            command = self.create_command(*args, **kwargs)
            logger.debug(f"Running command: {' '.join(command)}")
            subprocess.run(command, check=True)
            logger.info(f"Command completed successfully: {self.__class__.__name__}")
            if self.metadata_tracker:
                self._update_metadata(*args, **kwargs)
        except subprocess.CalledProcessError as e:
            logger.error(f"Error running GROMACS command: {e}")
            raise

    @abstractmethod
    def create_command(self, *args, **kwargs):
        pass

    @abstractmethod
    def metadata(self, *args, **kwargs):
        pass

    def _update_metadata(self, *args, **kwargs):
        metadata = self.metadata(*args, **kwargs)
        self.metadata_tracker.add_step(step_name="GROMACS", details=metadata)
