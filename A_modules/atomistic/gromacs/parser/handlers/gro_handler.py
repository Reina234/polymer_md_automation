from A_modules.atomistic.gromacs.parser.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
import re
from typing import List, Optional, Tuple


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

        super().__init__(store_top_line=True)  # Store the top line
        self.expected_columns = expected_columns
        self.num_atoms = 0
        self.atom_data = []
        self._box_dimensions = None

    @property
    def box_dimensions(self) -> List[float]:
        if self._box_dimensions is None:
            raise ValueError("Box dimensions have not been set.")
        return self._box_dimensions

    @box_dimensions.setter
    def box_dimensions(self, value: str):
        tokens = value.split()
        if len(tokens) != 3:
            raise ValueError(
                "Box dimensions must contain exactly three values (x, y, z)."
            )
        self._box_dimensions = [float(dim) for dim in tokens]

    def process(self, section):
        """
        Process a section of a .gro file to parse atom lines and extract box dimensions.
        """
        self.section = section
        lines = [line.strip() for line in section.lines if line.strip()]

        # Top line (title)
        self.top_line = lines.pop(0)

        # Number of atoms
        self.num_atoms = int(lines.pop(0))

        for line in lines:
            if len(self.atom_data) < self.num_atoms:
                normalized_tokens = self._normalize_atom_line(line)
                self.atom_data.append(normalized_tokens)
            else:
                self.box_dimensions = line

    def _normalize_atom_line(self, line: str) -> List:
        """
        Normalize a `.gro` atom line into separate tokens.
        Handles combined or separate residue number/name and atom name/index formats.
        """
        comment = None
        if ";" in line:
            content, comment = line.split(";", 1)
            line = content.strip()
            comment = comment.strip()

        tokens = line.split()

        # Handle minimum token count
        if len(tokens) < 5:
            raise ValueError(f"Invalid atom line format: {line}")

        # Determine if fields are combined or split
        residue_number, residue_name, tokens = self._parse_residue_field(tokens)
        atom_name, atom_index, tokens = self._parse_atom_field(tokens)

        # Extract coordinates (last 3 tokens)
        try:
            x, y, z = map(float, tokens[:3])
        except ValueError as e:
            raise ValueError(f"Invalid coordinate format in line: {line}. Details: {e}")

        return [residue_number, residue_name, atom_name, atom_index, x, y, z, comment]

    def _parse_residue_field(self, tokens: List[str]) -> Tuple[int, str, List[str]]:
        """
        Parse residue number and residue name, handling both combined and separate formats.

        Args:
            tokens (List[str]): List of tokens from the atom line.

        Returns:
            Tuple[int, str, List[str]]: Parsed residue number, residue name, and the remaining tokens.
        """
        # Extract the first token as the residue field
        residue_field = tokens.pop(0)

        # Regex to handle combined formats (e.g., 1UNL)
        match = re.match(r"^(\d+)([A-Za-z]+)$", residue_field)
        if match:
            residue_number = int(match.group(1))
            residue_name = match.group(2)
        else:
            # Handle separate formats (e.g., "1 UNL")
            try:
                residue_number = int(residue_field)
                residue_name = tokens.pop(0)  # Next token is the residue name
            except (ValueError, IndexError) as e:
                raise ValueError(
                    f"Invalid residue field or missing residue name: {residue_field}. Details: {e}"
                )
        return residue_number, residue_name, tokens

    def _parse_atom_field(self, tokens: List[str]) -> Tuple[str, int, List[str]]:
        """
        Parse atom name and atom index, handling combined or separate formats.
        """
        atom_name = tokens.pop(0)
        atom_index = int(tokens.pop(0))  # Assume atom index is always the next token
        return atom_name, atom_index, tokens

    @property
    def content(self) -> pd.DataFrame:
        """
        Return atom data as a Pandas DataFrame.
        """
        return pd.DataFrame(self.atom_data, columns=self.expected_columns)

    @content.setter
    def content(self, new_content: pd.DataFrame):
        if list(new_content.columns) != self.expected_columns:
            raise ValueError(
                "Columns of the DataFrame do not match the expected format."
            )
        self.atom_data = new_content.values.tolist()

    def _export_content(self) -> List[str]:
        """
        Export content to `.gro` format lines.
        """
        lines = []
        for row in self.atom_data:
            content = f"{row[0]:5}{row[1]:>5}{row[2]:>5}{row[3]:5}{row[4]:8.3f}{row[5]:8.3f}{row[6]:8.3f}"
            if row[7]:
                content += f" ; {row[7]}"
            lines.append(content)

        # Add box dimensions
        if self._box_dimensions:
            lines.append(" ".join(f"{dim:.6f}" for dim in self._box_dimensions))
        return lines
