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
    def extract_box_dimensions(content: List[str]) -> Optional[List[float]]:
        for line in content:
            if line.startswith("CRYST1"):
                try:
                    box_dimensions = [
                        float(line[6:15].strip()),
                        float(line[15:24].strip()),
                        float(line[24:33].strip()),
                    ]
                    logger.info(f"Extracted box dimensions: {box_dimensions}")

                    if any(
                        v in {0, None} or not isinstance(v, (int, float))
                        for v in box_dimensions
                    ):
                        logger.warning(
                            f"Invalid box dimensions in CRYST1 line: {line.strip()}"
                        )
                        return None

                    return box_dimensions
                except ValueError:
                    logger.warning(
                        f"Invalid box dimensions in CRYST1 line: {line.strip()}"
                    )
                    return None

        logger.warning("CRYST1 line not found; box dimensions are unavailable.")
        return None

    @staticmethod
    def extract_atoms(content: List[str]) -> List[Dict[str, float]]:
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
                    logger.warning(f"[!] Invalid atom data in line: {line.strip()}")
        logger.info(f"Extracted {len(atoms)} atoms from PDB file.")
        return atoms

    @staticmethod
    def get_atom_coordinates(atoms: List[Dict[str, float]]) -> List[List[float]]:
        return [[atom["x"], atom["y"], atom["z"]] for atom in atoms]

    @staticmethod
    def add_or_replace_box_dimensions(
        content: List[str], box_dimensions: List[float]
    ) -> List[str]:
        new_cryst1_line = f"CRYST1{box_dimensions[0]:9.3f}{box_dimensions[1]:9.3f}{box_dimensions[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
        modified_content = []
        cryst1_replaced = False

        for line in content:
            if line.startswith("CRYST1"):
                logger.info(
                    f"Replacing existing CRYST1 line with new dimensions: {box_dimensions}"
                )
                modified_content.append(new_cryst1_line)
                cryst1_replaced = True
            else:
                modified_content.append(line)

        if not cryst1_replaced:
            logger.info(f"Adding CRYST1 line with dimensions: {box_dimensions}")
            modified_content.insert(0, new_cryst1_line)  # Add CRYST1 line at the top

        return modified_content

    @staticmethod
    def add_comment(content: List[str], comment: str) -> List[str]:
        formatted_comment = f"{PDB_COMMENT_PREFIX} {comment.strip()}\n"
        logger.info(f"Adding comment to PDB content: {comment}")
        return [formatted_comment] + content
