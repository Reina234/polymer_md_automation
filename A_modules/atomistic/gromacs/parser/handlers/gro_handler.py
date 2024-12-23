from A_modules.atomistic.gromacs.parser.handlers.base_handler import (
    BaseHandler,
)
import pandas as pd
import re
from typing import List, Optional, Tuple
import logging

logger = logging.getLogger(__name__)


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
        """
        Returns the box dimensions as a list of three floats.
        """
        return self._box_dimensions

    @box_dimensions.setter
    def box_dimensions(self, box_dimensions: Optional[List[float]]):
        """
        Sets the box dimensions. Validates that the input is a list of three floats.
        Allows setting to None.
        """
        if box_dimensions is None:
            self._box_dimensions = None
            return

        if not isinstance(box_dimensions, list) or len(box_dimensions) != 3:
            raise ValueError("Box dimensions must be a list of three floats.")

        self._box_dimensions = [float(dim) for dim in box_dimensions]

    def process(self, section):
        """
        Process a .gro file section to parse atom lines and extract box dimensions.
        """
        self.section = section
        lines = [line.strip() for line in section.lines if line.strip()]

        # Top line (title)
        self.top_line = lines.pop(0)

        # Parse number of atoms
        if not lines:
            raise ValueError("Missing number of atoms line in the .gro file.")
        try:
            self.num_atoms = int(lines.pop(0))
        except ValueError:
            raise ValueError("Expected an integer for the number of atoms.")

        # Parse atom data and box dimensions
        for i, line in enumerate(lines):
            if len(self.atom_data) < self.num_atoms:
                normalized_tokens = self._normalize_atom_line(line)
                self.atom_data.append(normalized_tokens)
            else:
                # Parse box dimensions
                box_dims = self.parse_box_dimensions(line)
                if self.validate_box_dimensions(box_dims):
                    self.box_dimensions = box_dims
                else:
                    logger.warning(f"Invalid box dimensions: {line}")
                break

        # Validate that all atoms and box dimensions are parsed correctly
        if len(self.atom_data) != self.num_atoms:
            logger.error(
                f"Expected {self.num_atoms} atoms but parsed {len(self.atom_data)}."
            )
            raise ValueError("Mismatch between expected and parsed number of atoms.")

        if self.box_dimensions is None:
            logger.warning("No valid box dimensions found.")

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

    def _parse_box_dimensions(self, line: str) -> Optional[List[float]]:
        """
        Parse box dimensions from the last line of the .gro file.
        Returns None if the line does not contain valid box dimensions.
        """
        try:
            tokens = line.split()

            # Box dimensions must contain exactly 3 numeric tokens
            if len(tokens) == 3 and all(self._is_float(dim) for dim in tokens):
                return [float(dim) for dim in tokens]

            # If not valid, assume it's not a box dimensions line
            return None
        except ValueError as e:
            logger.warning(f"Error parsing box dimensions: {line}. Details: {e}")
            return None

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
        Export the parsed content to `.gro` format lines.
        """
        lines = []

        # Export atom data
        for row in self.atom_data:
            content = (
                f"{row[0]:5}{row[1]:>5}{row[2]:>5}{row[3]:5}"
                f"{row[4]:8.3f}{row[5]:8.3f}{row[6]:8.3f}"
            )
            if row[7]:
                content += f" ; {row[7]}"
            lines.append(content)

        # Export box dimensions if valid
        if self.validate_box_dimensions(self.box_dimensions):
            lines.append(" ".join(f"{dim:.6f}" for dim in self.box_dimensions))
        else:
            logger.error("Cannot export .gro file without valid box dimensions.")
            raise ValueError("Invalid or missing box dimensions during export.")

        return lines

    @staticmethod
    def _is_float(value: str) -> bool:
        """
        Check if a string can be converted to a float.
        """
        try:
            float(value)
            return True
        except ValueError:
            return False

    # NOTE: replace with the one in the utils functions !
    @staticmethod
    def parse_box_dimensions(line: str) -> Optional[List[float]]:
        """
        Parse a line containing box dimensions.
        Returns:
            List[float]: Parsed box dimensions if valid.
            None: If the line does not contain valid box dimensions.
        """
        try:
            tokens = line.split()
            if len(tokens) == 3 and all(GroHandler._is_float(dim) for dim in tokens):
                return [float(dim) for dim in tokens]
        except ValueError:
            pass
        return None

    @staticmethod
    def validate_box_dimensions(box_dimensions: Optional[List[float]]) -> bool:
        """
        Validate that box dimensions are a list of three positive floats.
        """
        if not box_dimensions or len(box_dimensions) != 3:
            return False
        return all(isinstance(dim, float) and dim > 0 for dim in box_dimensions)


# NOTE: formatting issues exist for the box dims part, so editconf is recommended instead :(
