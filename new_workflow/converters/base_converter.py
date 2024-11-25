from abc import ABC, abstractmethod

class FileConverter(ABC):
    """
    Abstract base class for file converters.
    """

    def __init__(self, metadata_manager: MetadataManager):
        self.metadata_manager = metadata_manager
        self.register_metadata()

    @abstractmethod
    def supported_formats(self) -> tuple:
        """
        Return (input_format, output_format) to indicate supported conversions.
        """
        pass

    @abstractmethod
    def convert(self, input_file: str, output_dir: str) -> str:
        """
        Convert the input file to the specified output format.
        """
        pass

    @abstractmethod
    def metadata(self) -> dict:
        """
        Provide metadata for this converter.
        """
        pass

    def register_metadata(self):
        """
        Register this converter's metadata with the MetadataManager.
        """
        data = self.metadata()
        self.metadata_manager.register_component(
            name=data["name"],
            version=data["version"],
            description=data["description"],
            options=data.get("options", {})
        )
