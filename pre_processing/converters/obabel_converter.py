from pre_processing.converters.base_file_converter import FileConverter
from pre_processing.metadata.metadata_manager import MetadataManager
import subprocess 
import os 

class OpenBabelConverter(FileConverter):
    """
    Converts molecule files using Open Babel.
    """

    def __init__(self, input_format: str, output_format: str, metadata_manager: MetadataManager):
        self.input_format = input_format
        self.output_format = output_format
        super().__init__(metadata_manager)

    def supported_formats(self) -> tuple:
        return (self.input_format, self.output_format)

    def convert(self, input_file: str, output_dir: str) -> str:
        output_file = os.path.join(
            output_dir, os.path.splitext(os.path.basename(input_file))[0] + f".{self.output_format}"
        )
        subprocess.run(
            [
                "obabel",
                input_file,
                "-O", output_file,
                "--gen3D"
            ],
            check=True
        )
        return output_file

    def metadata(self) -> dict:
        return {
            "name": "Open Babel Converter",
            "version": "3.1.1",
            "description": f"Converts files from {self.input_format} to {self.output_format} using Open Babel.",
            "options": {
                "3D coordinates": "Generated if not present",
                "tool": "Open Babel"
            }
        }
