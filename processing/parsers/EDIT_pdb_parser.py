from processing.parsers.base_parser import BaseParser
from typing import List, Optional
import logging
from config.constants import AVOGADROS_NUMBER, ANGSTROM_TO_CM
import os 


#NOTE: having a strict cubic definition could cause issues, e.g. if the pdb file has coords outside of it, is there another way of getting these box coords, like rdkit? 
logger = logging.getLogger(__name__)

class PDBParser(BaseParser):
    """
    A parser for handling PDB files, including extracting, validating,
    and modifying box information.
    """

    def __init__(self, file_path: str):
        super().__init__(file_path)
        self._content = self._read_file()
        self.box_dimensions = self._extract_box_dimensions()
        self.atoms = self._extract_atoms()
        self.num_molecules = self._count_molecules()

    def _read_file(self) -> List[str]:
        """
        Read the PDB file into memory.

        Returns:
            List[str]: The file content as a list of lines.
        """
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"PDB file not found: {self.file_path}")
        with open(self.file_path, "r") as file:
            return file.readlines()

    def _extract_box_dimensions(self) -> Optional[List[float]]:
        """
        Extract box dimensions from the CRYST1 line.

        Returns:
            Optional[List[float]]: Box dimensions [x, y, z] in Ångströms.

        Raises:
            ValueError: If any of the box dimensions are zero.
        """
        for line in self._content:
            if line.startswith("CRYST1"):
                try:
                    x = float(line[6:15].strip())
                    y = float(line[15:24].strip())
                    z = float(line[24:33].strip())

                    if x == 0 or y == 0 or z == 0:
                        raise ValueError(f"Box dimensions cannot be zero. Found: x={x}, y={y}, z={z}")

                    logger.info(f"[+] Extracted box dimensions: x={x}, y={y}, z={z}")
                    return [x, y, z]
                except ValueError as e:
                    raise ValueError(f"Invalid box dimensions in CRYST1 line: {e}")
        logger.warning("CRYST1 line not found; box dimensions are not available.")
        return None

    def _extract_atoms(self) -> List[dict]:
        """
        Extract atom information from the PDB file.

        Returns:
            List[dict]: A list of atom dictionaries containing atomic information.
        """
        atoms = []
        for line in self._content:
            if line.startswith(("ATOM", "HETATM")):
                atom = {
                    'serial': int(line[6:11].strip()),
                    'name': line[12:16].strip(),
                    'altLoc': line[16].strip(),
                    'resName': line[17:20].strip(),
                    'chainID': line[21].strip(),
                    'resSeq': int(line[22:26].strip()),
                    'iCode': line[26].strip(),
                    'x': float(line[30:38].strip()),
                    'y': float(line[38:46].strip()),
                    'z': float(line[46:54].strip()),
                    'element': line[76:78].strip(),
                }
                atoms.append(atom)
        logger.info(f"[+] Extracted {len(atoms)} atoms from PDB file.")
        return atoms

    def _count_molecules(self) -> int:
        """
        Count the number of unique molecules in the PDB file.

        Returns:
            int: Number of molecules.
        """
        # Assuming each molecule corresponds to a unique chain ID
        chain_ids = set(atom['chainID'] for atom in self.atoms if atom['chainID'])
        num_molecules = len(chain_ids) if chain_ids else 1  # Default to 1 if no chain IDs
        logger.info(f"[+] Number of molecules: {num_molecules}")
        return num_molecules

    def save(self, output_path: Optional[str] = None) -> str:
        """
        Save the PDB file to the specified output path.

        Args:
            output_path (Optional[str]): Path to save the modified file. Defaults to original path.

        Returns:
            str: Path to the saved file.
        """
        output_path = output_path or self.file_path
        with open(output_path, "w") as file:
            file.writelines(self._content)
        logger.info(f"[+] PDB file saved to {output_path}")
        return output_path

    def move(self, destination: str) -> None:
        """
        Move the file to a new location.

        Args:
            destination (str): Directory to move the file to.
        """
        if not os.path.exists(destination):
            os.makedirs(destination)
        new_path = os.path.join(destination, os.path.basename(self.file_path))
        os.rename(self.file_path, new_path)
        logger.info(f"[+] PDB file moved to {new_path}")
        self.file_path = new_path

    def has_box_information(self) -> bool:
        """
        Check if the PDB file includes box information (CRYST1 line).

        Returns:
            bool: True if the CRYST1 line is present, False otherwise.
        """
        return self.box_dimensions is not None

    def add_or_replace_box_information(self, box_size: List[float]) -> None:
        """
        Add or replace the box information in the PDB file.

        Args:
            box_size (List[float]): Dimensions of the box (x, y, z) in Ångströms.
        """
        box_line = f"CRYST1{box_size[0]:9.3f}{box_size[1]:9.3f}{box_size[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"

        # Check if CRYST1 exists and replace it
        for i, line in enumerate(self._content):
            if line.startswith("CRYST1"):
                self._content[i] = box_line
                self.box_dimensions = box_size
                logger.info("[+] Box information replaced in PDB file.")
                return

        # If no CRYST1 line, add it at the correct position
        insertion_index = next(
            (i for i, line in enumerate(self._content) if line.startswith(("ATOM", "HETATM", "MODEL", "END"))),
            len(self._content),
        )
        self._content.insert(insertion_index, box_line)
        self.box_dimensions = box_size
        logger.info("[+] Box information added to PDB file.")

    def calculate_density(self, molecular_weight: float) -> float:
        """
        Calculate the density of the PDB box.

        Args:
            molecular_weight (float): Molecular weight of the molecule in g/mol.

        Returns:
            float: Calculated density in g/cm³.

        Raises:
            ValueError: If box dimensions are not available.
        """
        if not self.box_dimensions:
            raise ValueError("Cannot calculate density without box dimensions.")

        # Total mass = number of molecules * molecular weight
        num_molecules = self.num_molecules

        total_mass_g = num_molecules * molecular_weight / AVOGADROS_NUMBER  # mass in grams

        # Calculate volume
        x, y, z = self.box_dimensions  # Box dimensions in Ångströms
        volume_cm3 = (x * ANGSTROM_TO_CM) * (y * ANGSTROM_TO_CM) * (z * ANGSTROM_TO_CM)  # Volume in cm³

        # Density = mass / volume
        density_g_cm3 = total_mass_g / volume_cm3

        logger.info(f"[+] Calculated density: {density_g_cm3:.6f} g/cm³")
        return density_g_cm3

    def add_comment(self, comment: str) -> None:
        """
        Add a comment to the PDB file.

        Args:
            comment (str): Comment to add to the top of the file.
        """
        formatted_comment = f"# {comment.strip()}\n"
        self._content.insert(0, formatted_comment)
        logger.info(f"[+] Comment added to PDB file: {comment}")
