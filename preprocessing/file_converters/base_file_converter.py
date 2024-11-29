from abc import ABC, abstractmethod
from typing import Optional
from preprocessing.metadata_tracker import MetadataTracker
import logging
import os
from preprocessing.utils import check_file_exists, check_file_type

# Set up logging (if needed for additional error tracking)
logging.basicConfig(level=logging.ERROR)

class BaseFileConverter(ABC):
    """
    Abstract base class for file converters.
    """
    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    @property
    @abstractmethod
    def input_file_type(self) -> str:
        """Return the expected input file type."""

    @property
    @abstractmethod
    def output_file_type(self) -> str:
        """Return the expected output file type."""

    @abstractmethod
    def convert(self, input_file_path: str, output_dir: Optional[str]= None) -> str:
        pass

    @property 
    @abstractmethod 
    def converter_method(self) -> str:
        """Return the method"""

   
    def get_output_file_path(self, input_file_path: str, output_dir: Optional[str] = None) -> str:
        """
        Generates the output file name based on the input file name,
        replacing its extension with the expected output file type
        and placing it in the specified output directory.

        Args:
            input_file_path (str): The path to the input file.
            output_dir (str): The directory for the output file.

        Returns:
            str: The full path of the output file.
        """

        if output_dir is None:
            output_dir = os.path.dirname(input_file_path)
        check_file_exists(input_file_path)
        check_file_type(input_file_path, self.input_file_type)

        base_name = os.path.splitext(os.path.basename(input_file_path))[0]
        output_file_name = f"{base_name}.{self.output_file_type}"
        output_file_path = os.path.join(output_dir, output_file_name)

        return output_file_path

    def _log_conversion(self, input_file_path: str, output_file_path: str) -> None:
        """
        Log the conversion of an input file to an output file.
        """
        logging.info(f"Converted {input_file_path} to {output_file_path} using {self.converter_method}.")


    def metadata(self, input_file_path: str, output_file_path: str, additional_notes:Optional[str] = None) -> dict:
        """
        Return metadata describing the conversion tool and its configuration.
        """
        return {"program(s) used": self.converter_method,
                "details": f".{self.input_file_type} -> .{self.output_file_type}",
                "action(s)": f"converted input file at {input_file_path}, saved to {output_file_path}",
                "additional note(s)": additional_notes}

    def _update_metadata(self, input_file_path: str, output_file_path: str, additional_notes: Optional[str] = None) -> None:
        """
        Add metadata to the metadata tracker.
        """
        metadata = self.metadata(input_file_path, output_file_path, additional_notes)
        self.metadata_tracker.add_step(step_name="Conversion", details=metadata)

#NOTE: should I have an if here? 