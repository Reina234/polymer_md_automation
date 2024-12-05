from ZZZ_parser.handlers.base_handler import BaseHandler
from typing import List


class IncludesHandler(BaseHandler):
    def __init__(self):
        super().__init__(store_top_line=False)  # No static top line for #include
        self.include_line = None  # Stores the #include line

    def process_line(self, line: str):
        """
        Processes the #include line.
        """
        self.include_line = line.strip()

    @property
    def content(self) -> str:
        """
        Returns the #include line (with in-line comments).
        """
        return self.include_line

    @content.setter
    def content(self, new_content: str):
        """
        Sets the #include line.
        """
        self.include_line = new_content

    def _export_content(self) -> List[str]:
        """
        Exports the #include line.
        """
        return [self.include_line]
