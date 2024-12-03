from typing import List, Optional
import os
from abc import ABC
from A_modules.shared.utils import check_file_exists, get_file_contents, check_file_type
from A_modules.shared.metadata_tracker import MetadataTracker
from abc import abstractmethod, ABC


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
