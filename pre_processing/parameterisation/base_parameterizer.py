from abc import ABC, abstractmethod
from pre_processing.converters.conversion_service import ConversionService
from pre_processing.metadata.metadata_manager import MetadataManager  

class BaseParameterizer(ABC):
    """
    Abstract base class for parameterization methods.
    """

    def __init__(self, method_name: str, conversion_service: ConversionService, metadata_manager: MetadataManager):
        self.method_name = method_name
        self.conversion_service = conversion_service
        self.metadata_manager = metadata_manager
        self.register_metadata()

    @abstractmethod
    def parameterize(self, input_file: str, output_dir: str) -> str:
        pass

    @abstractmethod
    def metadata(self) -> dict:
        """
        Provide metadata for this parameterizer.
        """
        pass

    def register_metadata(self):
        """
        Register this parameterizer's metadata with the MetadataManager.
        """
        data = self.metadata()
        self.metadata_manager.register_component(
            name=data["name"],
            version=data["version"],
            description=data["description"],
            options=data.get("options", {})
        )
