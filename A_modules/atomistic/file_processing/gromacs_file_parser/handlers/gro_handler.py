from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
from typing import List


class GroHandler(BaseHandler):
    re_pattern = None
    construct_name = "gro"
    suppress = None

    def __init__(
        self,
        expected_columns: List[str] = [
            "Residue Number",
            "Residue Name",
            "Atom Name",
            "Atom Index",
            "X",
            "Y",
            "Z",
            "In-line comments",
        ],
    ):
        self.expected_columns = expected_columns
        super().__init__(store_top_line=True)  # Store the top line
        self.num_atoms = 0  # Number of atoms
        self.atom_data = []  # To store atom rows
        self._box_dimensions = None  # To store box dimensions as a string

    @property
    def box_dimensions(self) -> List[float]:
        """
        Returns the box dimensions as a list of floats.
        """
        if self._box_dimensions is None:
            raise ValueError("Box dimensions have not been set.")
        return self._box_dimensions

    @box_dimensions.setter
    def box_dimensions(self, value: str):
        """
        Validates and sets the box dimensions from a space-separated string.
        The input must be in the format: "x y z".

        :param value: A space-separated string representing the box dimensions.
        :type value: str
        """
        tokens = value.split()
        if len(tokens) != 3:
            raise ValueError(
                "Box dimensions must contain exactly three values (x, y, z)."
            )

        try:
            self._box_dimensions = [float(dim) for dim in tokens]
        except ValueError:
            raise ValueError("Box dimensions must be valid floating-point numbers.")

    def process(self, section):
        """
        Extends the process method to handle the box dimensions (bottom line).
        """
        self.section = section

        # Filter out empty lines
        lines = [line.strip() for line in section.lines if line.strip()]

        # Handle the top line
        if self.store_top_line and lines:
            self.top_line = lines.pop(0)

        # Handle number of atoms
        if lines:
            self.num_atoms = int(
                lines.pop(0)
            )  # First remaining line is number of atoms

        # Handle atom data and box dimensions
        for line in lines:
            if len(self.atom_data) < self.num_atoms:
                self.atom_data.append(self._parse_atom_line(line))
            else:
                # Last line after atom data is box dimensions
                self.box_dimensions = line

    def _parse_atom_line(self, line: str) -> List:
        """
        Parses a single line of atom data in the .gro file, handling variable spacing and in-line comments.

        :param line: A line of atom data from the .gro file.
        :type line: str
        :return: Parsed atom data as a list.
        :rtype: List
        """
        if ";" in line:
            content, comment = line.split(";", 1)
            line = content.strip()
            comment = comment.strip()
        else:
            comment = None

        # Split the line by whitespace
        tokens = line.split()

        if len(tokens) < 7:
            raise ValueError(f"Invalid atom line format: {line}")

        # Parse the tokens
        residue_number = int(tokens[0])  # Residue number
        residue_name = tokens[1]  # Residue name
        atom_name = tokens[2]  # Atom name
        atom_index = int(tokens[3])  # Atom index
        x = float(tokens[4])  # X coordinate
        y = float(tokens[5])  # Y coordinate
        z = float(tokens[6])  # Z coordinate

        return [residue_number, residue_name, atom_name, atom_index, x, y, z, comment]

    @property
    def content(self) -> pd.DataFrame:
        """
        Returns the atom data as a DataFrame.
        """
        columns = self.expected_columns
        return pd.DataFrame(self.atom_data, columns=columns)

    @content.setter
    def content(self, new_content: pd.DataFrame):
        """
        Updates atom data using a DataFrame.
        """
        expected_columns = self.expected_columns
        if list(new_content.columns) != expected_columns:
            raise ValueError(
                "Columns of the DataFrame do not match the expected atom data format."
            )
        self.atom_data = new_content.values.tolist()

    def _export_content(self) -> List[str]:
        """
        Exports the .gro file content as a list of lines, including box dimensions.
        """
        lines = []

        # Atom data
        for row in self.atom_data:
            content = f"{row[0]:5}{row[1]:>5}{row[2]:>5}{row[3]:5}{row[4]:8.3f}{row[5]:8.3f}{row[6]:8.3f}"
            inline_comment = row[7] if len(row) > 7 else None
            if inline_comment:
                lines.append(f"{content} ; {inline_comment}")
            else:
                lines.append(content)

        # Add box dimensions at the end
        if self._box_dimensions:
            lines.append(" ".join(f"{dim:.6f}" for dim in self._box_dimensions))

        return lines
