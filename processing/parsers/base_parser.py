# File: parsers/base_parser.py
import os
from typing import List
from abc import ABC, abstractmethod
from typing import Optional

class BaseParser(ABC):
    """
    Abstract base class for file parsers.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path

    def _read_file(self) -> List[str]:
        """
        Read the PDB file into memory.

        Returns:
            List[str]: The file content as a list of lines.
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"PDB file not found: {self.file_path}")
        with open(self.file_path, "r") as file:
            return file.readlines()
        
    @abstractmethod
    def save(self, output_path: Optional[str] = None):
        pass

    @abstractmethod
    def move(self, destination: str):
        pass

    @abstractmethod
    def add_comment(self, comment: str):
        pass
