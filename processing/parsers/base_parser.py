# File: parsers/base_parser.py

from abc import ABC, abstractmethod
from typing import Optional

class BaseParser(ABC):
    """
    Abstract base class for file parsers.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path

    @abstractmethod
    def _read_file(self):
        pass

    @abstractmethod
    def save(self, output_path: Optional[str] = None):
        pass

    @abstractmethod
    def move(self, destination: str):
        pass

    @abstractmethod
    def add_comment(self, comment: str):
        pass
