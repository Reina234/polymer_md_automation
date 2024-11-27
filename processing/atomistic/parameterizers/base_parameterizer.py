# File: processing/common/base_parameterizer.py

from abc import ABC, abstractmethod
from processing.metadata_tracker import MetadataTracker
from processing.metadata_tracker import MetadataTracker


class BaseParameterizer(ABC):
    """
    Abstract base class for parameterizers.
    """
    BASE_TEMPORARY_DIRECTORY = "temp/atomistic"
    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    @abstractmethod
    def parameterize(self, input_dir: str, output_dir: str) -> dict:
        """
        Abstract method to parameterize input data.
        Args:
            input_dir (str): Directory containing the input files.
            output_dir (str): Directory to save output files.
        Returns:
            dict: Paths to the generated files.
        """
        pass

    @abstractmethod
    def metadata(self) -> dict:
        """
        Return metadata about the parameterizer.
        """
        pass
