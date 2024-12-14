# File: processing/common/base_parameterizer.py

from abc import ABC, abstractmethod
from A_modules.shared.metadata_tracker import MetadataTracker
from typing import Optional


class BaseParameterizer(ABC):
    """
    Abstract base class for parameterizers.
    """

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
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

    def _update_metadata(
        self,
        input_file_path: str,
        output_dir: str,
        additional_notes: Optional[str] = None,
    ) -> None:
        """
        Add metadata to the metadata tracker.
        """
        metadata = self.metadata(input_file_path, output_dir, additional_notes)
        self.metadata_tracker.add_step(step_name="Parameterization", details=metadata)

    @abstractmethod
    def metadata(
        self,
        input_file_path: str,
        output_dir: str,
        additional_notes: Optional[str] = None,
    ) -> dict:
        """
        Return metadata about the parameterizer.
        """
        pass
