class ConversionService:
    """
    Manages file conversions and delegates tasks to registered converters.
    """

    def __init__(self):
        self._registry = []

    def register_converter(self, converter):
        self._registry.append(converter)

    def convert(self, input_file: str, output_format: str, output_dir: str) -> str:
        input_format = input_file.split('.')[-1]
        for converter in self._registry:
            if converter.supported_formats() == (input_format, output_format):
                return converter.convert(input_file, output_dir)
        raise ValueError(f"No converter found for {input_format} to {output_format}")
