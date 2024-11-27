import os
import logging
from typing import List, Optional, Dict

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class PDBParser:
    """
    A class to parse PDB files and extract relevant information.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.content = self._read_file()
        self.box_dimensions = self._extract_box_dimensions()  # Box dimensions in nm
        self.atoms = self._extract_atoms()  # Atom data

    def _read_file(self) -> List[str]:
        """
        Read the PDB file and return its content as a list of lines.
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"PDB file not found: {self.file_path}")

        with open(self.file_path, "r") as file:
            lines = file.readlines()

        if not lines:
            raise ValueError(f"PDB file is empty: {self.file_path}")

        logger.info(f"[+] PDB file read successfully: {self.file_path}")
        return lines

    def _extract_box_dimensions(self) -> Optional[List[float]]:
        """
        Extract box dimensions from the CRYST1 line in the PDB file.

        Returns:
            Optional[List[float]]: Box dimensions [x, y, z] in nm.
        """
        for line in self.content:
            if line.startswith("CRYST1"):
                try:
                    x = float(line[6:15].strip())
                    y = float(line[15:24].strip())
                    z = float(line[24:33].strip())
                    logger.info(f"[+] Extracted box dimensions: {x}, {y}, {z}")
                    return [x, y, z]
                except ValueError:
                    logger.warning(f"[!] Invalid box dimensions in CRYST1 line: {line.strip()}")
        logger.warning("[!] CRYST1 line not found; box dimensions are unavailable.")
        return None

    def _extract_atoms(self) -> List[Dict]:
        """
        Extract atom information from the PDB file.

        Returns:
            List[Dict]: List of atoms with their x, y, z coordinates.
        """
        atoms = []
        for line in self.content:
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

    def get_atom_coordinates(self) -> List[List[float]]:
        """
        Get atom coordinates as a list of [x, y, z].

        Returns:
            List[List[float]]: List of atom coordinates.
        """
        return [[atom["x"], atom["y"], atom["z"]] for atom in self.atoms]

    def save(self, output_path: Optional[str] = None) -> str:
        """
        Save the PDB file to the specified output path.

        Args:
            output_path (Optional[str]): Path to save the modified file. Defaults to the original file path.

        Returns:
            str: Path to the saved file.
        """
        output_path = output_path or self.file_path
        with open(output_path, "w") as file:
            file.writelines(self.content)
        logger.info(f"[+] PDB file saved to {output_path}")
        return output_path
