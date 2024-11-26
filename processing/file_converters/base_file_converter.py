from abc import ABC, abstractmethod
from typing import Tuple


class BaseFileConverter(ABC):
    """
    Abstract base class for file converters.
    """

    @abstractmethod
    def supported_formats(self) -> Tuple[str, str]:
        """
        Return a tuple (input_format, output_format) to indicate the supported conversion.
        Example: ("pdb", "mol2")
        """
        pass

    @abstractmethod
    def convert(self, input_file: str, output_dir: str) -> str:
        """
        Convert the input file to the desired format.
        Parameters:
        - input_file: Path to the input file.
        - output_dir: Directory to save the converted file.
        Returns:
        - Path to the converted output file.
        """
        pass

    @abstractmethod
    def metadata(self) -> dict:
        """
        Return metadata describing the conversion tool and its configuration.
        """
        pass
