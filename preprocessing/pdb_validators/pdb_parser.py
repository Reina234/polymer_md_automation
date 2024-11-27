import os
import logging
from typing import List, Optional, Dict
from preprocessing.utils import check_file_exists, check_file_type

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


import os
import logging
from typing import List, Optional, Dict

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class PDBParser:
    """
    A stateless class to parse PDB files and extract relevant information.
    """

    @staticmethod
    def read_file(input_file_path: str) -> List[str]:
        """
        Read the PDB file and return its content as a list of lines.

        Args:
            input_file_path (str): Path to the PDB file.

        Returns:
            List[str]: Lines from the PDB file.
        """
        check_file_exists(input_file_path)
        content = check_file_exists(input_file_path, "pdb")
        return content

    @staticmethod
    def extract_box_dimensions(content: List[str]) -> Optional[List[float]]:
        """
        Extract box dimensions from the CRYST1 line in the PDB file.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            Optional[List[float]]: Box dimensions [x, y, z] in nm, or None if not found.
        """
        for line in content:
            if line.startswith("CRYST1"):
                try:
                    x = float(line[6:15].strip())
                    y = float(line[15:24].strip())
                    z = float(line[24:33].strip())
                    logger.info(f"[+] Extracted box dimensions: {x}, {y}, {z}")
                    return [x, y, z]
                except ValueError:
                    logger.warning(
                        f"[!] Invalid box dimensions in CRYST1 line: {line.strip()}"
                    )
        logger.warning("[!] CRYST1 line not found; box dimensions are unavailable.")
        return None

    @staticmethod
    def extract_atoms(content: List[str]) -> List[Dict[str, float]]:
        """
        Extract atom information from the PDB file.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            List[Dict[str, float]]: List of atoms with their x, y, z coordinates.
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
                    logger.warning(f"[!] Invalid atom data in line: {line.strip()}")
        logger.info(f"[+] Extracted {len(atoms)} atoms from PDB file.")
        return atoms

    @staticmethod
    def get_atom_coordinates(atoms: List[Dict[str, float]]) -> List[List[float]]:
        """
        Get atom coordinates as a list of [x, y, z].

        Args:
            atoms (List[Dict[str, float]]): List of atom data.

        Returns:
            List[List[float]]: List of atom coordinates.
        """
        return [[atom["x"], atom["y"], atom["z"]] for atom in atoms]

    @staticmethod
    def save(output_file_path: str, content: List[str]) -> str:
        """
        Save the PDB file to the specified file path.

        Args:
            output_file_path (str): Path to save the file.
            content (List[str]): Lines to write to the file.

        Returns:
            str: Path to the saved file.
        """
        with open(output_file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] PDB file saved to {output_file_path}")
        return output_file_path

    @staticmethod
    def add_or_replace_box_dimensions(
        content: List[str], box_dimensions: List[float]
    ) -> List[str]:
        """
        Add or replace the CRYST1 line in the PDB content with the given box dimensions.

        Args:
            content (List[str]): The content of the PDB file as a list of lines.
            box_dimensions (List[float]): The box dimensions [x, y, z] in nm.

        Returns:
            List[str]: The modified content with the updated CRYST1 line.
        """
        new_cryst1_line = f"CRYST1{box_dimensions[0]:9.3f}{box_dimensions[1]:9.3f}{box_dimensions[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
        modified_content = []
        cryst1_replaced = False

        for line in content:
            if line.startswith("CRYST1"):
                logger.info(
                    f"[+] Replacing existing CRYST1 line with new dimensions: {box_dimensions}"
                )
                modified_content.append(new_cryst1_line)
                cryst1_replaced = True
            else:
                modified_content.append(line)

        if not cryst1_replaced:
            logger.info(f"[+] Adding CRYST1 line with dimensions: {box_dimensions}")
            modified_content.insert(0, new_cryst1_line)  # Add CRYST1 line at the top

        return modified_content

    @staticmethod
    def add_comment(content: List[str], comment: str) -> List[str]:
        """
        Add a comment to the top of the PDB content.

        Args:
            content (List[str]): The content of the PDB file as a list of lines.
            comment (str): The comment to add.

        Returns:
            List[str]: The modified content with the comment added at the top.
        """
        formatted_comment = f"# {comment.strip()}\n"
        logger.info(f"[+] Adding comment to PDB file: {comment}")
        return [formatted_comment] + content
