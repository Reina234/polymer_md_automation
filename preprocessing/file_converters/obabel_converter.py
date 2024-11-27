import os
import subprocess
from preprocessing.file_converters.base_file_converter import BaseFileConverter
from preprocessing.metadata_tracker import MetadataTracker
from typing import Optional

class OpenBabelConverter(BaseFileConverter):
    """
    Converts molecule files using Open Babel.
    """

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    @property
    def input_file_type(self) -> str:
        return "pdb"
    
    @property
    def output_file_type(self) -> str:
        return "mol2"
    
    @property
    def converter_method(self) -> str:
        return "Open Babel"

    def convert(self, input_file_path: str, output_dir: Optional[str] = None, additional_notes: Optional[str] = None) -> str:
        """
        Convert the input file to the specified output format.
        """

        output_file_path = self.get_output_file_path(input_file_path, output_dir)

        try:
            subprocess.run(
                [
                    "obabel",
                    input_file_path,
                    "-O", output_file_path,
                    "--gen3D",
                      "--h"  # add 3D coords via --h
                ],
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Open Babel failed to convert {input_file_path} to {self.output_file_type}. Error: {e}")

        if self.metadata_tracker is not None:
            self._update_metadata(input_file_path, output_file_path, additional_notes)

        self._log_conversion(input_file_path, output_file_path)

        return output_file_path

