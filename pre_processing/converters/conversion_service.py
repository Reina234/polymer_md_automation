from pre_processing.converters.base_file_converter import BaseFileConverter
import os 

class ConversionService:
    """
    Manages file conversions and delegates tasks to registered converters.
    """

    def __init__(self):
        self._registry = []

    def register_converter(self, converter: BaseFileConverter):
        """
        Register a new file converter.
        """
        self._registry.append(converter)

    def convert(self, input_file: str, output_format: str, output_dir: str) -> str:
        """
        Find an appropriate converter and perform the conversion.
        """
        input_format = os.path.splitext(input_file)[1][1:]  # Extract input file extension
        for converter in self._registry:
            if converter.supported_formats() == (input_format, output_format):
                return converter.convert(input_file, output_dir)
        raise ValueError(f"No converter found for {input_format} to {output_format}")

    def list_supported_conversions(self) -> list:
        """
        List all supported conversions.
        """
        return [converter.supported_formats() for converter in self._registry]
