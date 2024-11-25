from abc import ABC, abstractmethod
from typing import Optional
from solvent.solvent import Solvent
class PDBProvider(ABC):
    """
    Abstract base class for generating or loading PDB files from various sources.
    """

    def __init__(self, default_note: str, solvent: Solvent, additional_notes: Optional[str] = None):
        """
        :param default_note: Default note describing the PDB source.
        :param additional_notes: Optional user-provided notes.
        """
        self.default_note = default_note
        self._additional_notes = additional_notes  # Encapsulated attribute
        self.solvent = solvent
        self._molecule_name = solvent.name

    @property
    def additional_notes(self) -> Optional[str]:
        """
        Getter for additional notes.
        """
        return self._additional_notes

    @additional_notes.setter
    def additional_notes(self, notes: str) -> None:
        """
        Setter for additional notes.
        """
        self._additional_notes = notes

    @property
    def metadata(self) -> str:
        """
        Retrieve the metadata for the PDB file, including default and additional notes.
        """
        if self.additional_notes:
            return f"{self.default_note}\nAdditional notes: {self.additional_notes}"
        return self.default_note

    @abstractmethod
    def get_pdb(self, output_dir: str, identifier: str = "") -> str:
        pass

    @abstractmethod
    def get_source(self) -> str:
        pass
