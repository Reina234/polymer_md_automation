from typing import List, Optional
import os
from abc import ABC
from A_modules.shared.utils import check_file_exists, get_file_contents, check_file_type
from A_modules.shared.metadata_tracker import MetadataTracker
from abc import abstractmethod, ABC
import logging

logger = logging.getLogger(__name__)


class BaseFileParser(ABC):
    """
    A base class to parse files
    """

    @abstractmethod
    def read_file(file_path: str, file_type: str) -> List[str]:
        pass

    @abstractmethod
    def add_comment(self, content: List[str], comments: List[str]) -> List[str]:
        pass


class FileParser:
    """
    A base class for parsing different file formats.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.content = self.read_file()

    def read_file(self) -> List[str]:
        """
        Reads the content of the file and returns it as a list of strings.
        """
        if not os.path.exists(self.file_path):
            logger.error(f"File not found: {self.file_path}")
            raise FileNotFoundError(f"File not found: {self.file_path}")

        with open(self.file_path, "r") as file:
            content = file.readlines()

        logger.info(f"Read {len(content)} lines from {self.file_path}")
        return content

    def parse(self):
        pass

    def validate(self):
        pass
