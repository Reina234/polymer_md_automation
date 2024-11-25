# core/pdb_provider.py
from abc import ABC, abstractmethod
from typing import Optional


class PDBProvider(ABC):
    """
    Abstract base class for providing PDB files.
    """

    @abstractmethod
    def get_pdb(self, output_file: str) -> None:
        """
        Generate or load a PDB file and save it to the specified path.
        :param output_file: Path to save the PDB file.
        """
        pass
