from abc import ABC, abstractmethod
from typing import Dict


class BaseParameterizer(ABC):
    """
    Abstract base class for parameterization methods.
    Defines a consistent interface for parameterizing molecular structures.
    """

    def __init__(self, metadata_tracker):
        """
        Initialize the parameterizer with a shared metadata tracker.
        """
        self.metadata_tracker = metadata_tracker

    @abstractmethod
    def parameterize(self, input_file: str, output_dir: str) -> Dict[str, str]:
        """
        Abstract method to parameterize a molecule.
        Parameters:
        - input_file: Path to the input structure file (e.g., MOL2).
        - output_dir: Directory to save parameterization outputs.
        Returns:
        - Dictionary with paths to generated parameterized files (e.g., GRO, TOP).
        """
        pass

    @abstractmethod
    def metadata(self) -> Dict:
        """
        Return metadata describing the parameterization tool and its configuration.
        """
        pass
