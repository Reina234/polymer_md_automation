import os
import subprocess
from processing.file_converters.base_file_converter import BaseFileConverter
from processing.metadata_tracker import MetadataTracker
from typing import Tuple

class OpenBabelConverter(BaseFileConverter):
    """
    Converts molecule files using Open Babel.
    """

    def __init__(self, input_format: str, output_format: str, metadata_tracker: MetadataTracker):
        self.input_format = input_format
        self.output_format = output_format
        self.metadata_tracker = metadata_tracker

    def supported_formats(self) -> Tuple[str, str]:
        """
        Return the input and output formats supported by this converter.
        """
        return self.input_format, self.output_format

    def convert(self, input_file: str, output_dir: str) -> str:
        """
        Convert the input file to the specified output format.
        """
        if not input_file.endswith(f".{self.input_format}"):
            raise ValueError(f"Expected a {self.input_format} file. Got: {input_file}")

        # Construct the output file path
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = os.path.join(output_dir, f"{base_name}.{self.output_format}")

        # Run Open Babel conversion
        try:
            subprocess.run(
                [
                    "obabel",
                    input_file,
                    "-O", output_file,
                    "--gen3D",
                      "--h"  # Add 3D coordinates if not present
                ],
                check=True
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"Open Babel failed to convert {input_file} to {self.output_format}. Error: {e}")

        # Update metadata
        self.metadata_tracker.add_step("Open Babel",{
                "input_file": input_file,
                "output_format": self.output_format,
                "output_file": output_file,
                "version": "3.1.1",
                "description": f"Converts {self.input_format} files to {self.output_format} files.",
            }
        )
        


        return output_file

    def metadata(self) -> dict:
        """
        Return metadata for the Open Babel conversion process.
        """
        return {
            "name": "Open Babel Converter",
            "version": "3.1.1",
            "description": f"Converts files from {self.input_format} to {self.output_format}.",
            "options": {
                "input_format": self.input_format,
                "output_format": self.output_format,
            }
        }
