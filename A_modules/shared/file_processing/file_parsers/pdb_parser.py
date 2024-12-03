import os
import logging
from typing import List, Optional, Dict
from A_modules.shared.file_processing.file_parsers.base_file_parser import (
    BaseFileParser,
)
from A_modules.shared.utils import check_file_exists, check_file_type, get_file_contents

logger = logging.getLogger(__name__)

PDB_FILE_EXTENSION = "pdb"
PDB_COMMENT_PREFIX = "#"


class PDBParser(BaseFileParser):
    """
    A stateless class to parse PDB files and extract relevant information.
    """

    @staticmethod
    def read_file(file_path: str, file_type=PDB_COMMENT_PREFIX) -> List[str]:
        file_path = check_file_exists(file_path)
        check_file_type(file_path, file_type)
        content = get_file_contents(file_path)
        return content

    def extract_box_dimensions(content: List[str]) -> Optional[List[float]]:
        """
        Extracts box dimensions from the CRYST1 line in the PDB file.

        :param content: Lines from the PDB file.
        :return: Box dimensions [x, y, z] or None if not found.
        """
        for line in content:
            if line.startswith("CRYST1"):
                try:
                    box_dimensions = [
                        float(line[6:15].strip()),
                        float(line[15:24].strip()),
                        float(line[24:33].strip()),
                    ]
                    logger.info(f"Extracted box dimensions: {box_dimensions}")
                    return box_dimensions
                except ValueError:
                    logger.warning(
                        f"Invalid box dimensions in CRYST1 line: {line.strip()}"
                    )
                    return None
        logger.warning("CRYST1 line not found; box dimensions unavailable.")
        return None

    def extract_atoms(content: List[str]) -> List[Dict[str, float]]:
        """
        Extracts atom information from the PDB file.

        :param content: Lines from the PDB file.
        :return: List of atoms with x, y, z coordinates.
        """
        atoms = []
        for line in content:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom = {
                        "x": float(line[30:38].strip()),
                        "y": float(line[38:46].strip()),
                        "z": float(line[46:54].strip()),
                    }
                    atoms.append(atom)
                except ValueError:
                    logger.warning(f"Invalid atom data in line: {line.strip()}")
        logger.info(f"Extracted {len(atoms)} atoms from PDB file.")
        return atoms

    @staticmethod
    def get_atom_coordinates(atoms: List[Dict[str, float]]) -> List[List[float]]:
        """
        Extracts atom coordinates as a list of [x, y, z].

        :param atoms: List of atoms with x, y, z coordinates.
        :return: List of atom coordinates.
        """
        if not atoms:
            logger.warning("No atoms found to extract coordinates.")
            return []

        # NOTE: need to debug this further
        if all(
            isinstance(atom, dict) and "x" in atom and "y" in atom and "z" in atom
            for atom in atoms
        ):
            return [[atom["x"], atom["y"], atom["z"]] for atom in atoms]

        logger.error("Atom data is not in the expected format.")
        return []

    def add_or_replace_box_dimensions(
        content: List[str], box_dimensions: List[float]
    ) -> List[str]:
        """
        Adds or replaces the CRYST1 line with the provided box dimensions.

        :param content: Lines from the PDB file.
        :param box_dimensions: Box dimensions [x, y, z].
        :return: Modified content with updated CRYST1 line.
        """
        new_cryst1_line = f"CRYST1{box_dimensions[0]:9.3f}{box_dimensions[1]:9.3f}{box_dimensions[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
        modified_content = []
        cryst1_replaced = False

        for line in content:
            if line.startswith("CRYST1"):
                logger.info(
                    f"Replacing CRYST1 line with new dimensions: {box_dimensions}"
                )
                modified_content.append(new_cryst1_line)
                cryst1_replaced = True
            else:
                modified_content.append(line)

        if not cryst1_replaced:
            logger.info(f"Adding CRYST1 line with dimensions: {box_dimensions}")
            modified_content.insert(0, new_cryst1_line)

        return modified_content

    def add_comment(content: List[str], comment: str) -> List[str]:
        """
        Adds a comment to the top of the PDB content.

        :param content: Lines from the PDB file.
        :param comment: The comment to add.
        :return: Modified content with comment added.
        """
        formatted_comment = f"{PDB_COMMENT_PREFIX} {comment.strip()}\n"
        logger.info(f"Adding comment to PDB file: {comment}")
        return [formatted_comment] + content
