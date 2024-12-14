from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
from typing import List


class DefaultHandler(BaseHandler):
    re_pattern = None
    construct_name = "default"
    suppress = None

    def __init__(self):
        super().__init__(store_top_line=False)  # Does not store a top line
        self.data = []  # Stores raw data rows
        self.expected_headers = []  # Can be set dynamically

    def process_line(self, line: str):
        """
        Processes a data line, splitting in-line comments if present.
        """
        if ";" in line:
            content, comment = line.split(";", 1)
            self.data.append(content.strip().split() + [comment.strip()])
        else:
            self.data.append(
                line.split() + [None]
            )  # Add None for missing in-line comment

    @property
    def content(self) -> pd.DataFrame:
        """
        Returns the data block as a DataFrame.
        """
        if not self.expected_headers:
            raise ValueError("Headers must be defined.")
        columns = self.expected_headers + ["In-Line Comments"]
        return pd.DataFrame(self.data, columns=columns)

    @content.setter
    def content(self, new_content: pd.DataFrame):
        """
        Sets the content using a DataFrame.
        """
        if list(new_content.columns) != self.expected_headers + ["In-Line Comments"]:
            raise ValueError("Headers do not match the expected format.")
        self.data = new_content.values.tolist()

    def _export_content(self) -> List[str]:
        """
        Exports the data block as lines, appending in-line comments where applicable.
        """
        lines = []
        for row in self.data:
            content = " ".join(row[:-1])  # All columns except the last
            inline_comment = row[-1]
            if inline_comment:
                lines.append(f"{content} ; {inline_comment}")
            else:
                lines.append(content)
        return lines
