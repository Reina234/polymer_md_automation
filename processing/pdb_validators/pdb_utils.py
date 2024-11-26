# File: utils/pdb_utils.py
from config.constants import AVOGADROS_NUMBER, CM_3_TO_NM_3_CONVERSION
from typing import List
# File: utils/pdb_utils.py

from typing import List
import math


def calculate_box_size(molecular_weight: float, density: float, target_molecules: int = 1000) -> List[float]:
    """
    Calculate the box size required for a given solvent based on its molecular weight and density.
    
    Args:
        molecular_weight (float): Molecular weight of the solvent in g/mol.
        density (float): Density of the solvent in g/cm³.
        target_molecules (int): Desired number of molecules in the box.

    Returns:
        List[float]: Box dimensions [x, y, z] in nm.
    """

    
      # Molecules per mol
    mass_per_molecule = molecular_weight / AVOGADROS_NUMBER  # g/molecule

    total_mass = target_molecules * mass_per_molecule  # g
    volume_cm3 = total_mass / density  # cm³

    volume_nm3 = volume_cm3 * CM_3_TO_NM_3_CONVERSION # Convert cm³ to nm³
    side_length = math.cbrt(volume_nm3)  # Cube root to find box side length in nm

    return [side_length, side_length, side_length]


def add_box_information(file_path: str, box_size: List[float]) -> None:
    """
    Add box information (CRYST1 line) to a PDB file in the correct position.
    
    Args:
        file_path (str): Path to the PDB file.
        box_size (List[float]): Dimensions of the box (x, y, z).
    """
    box_line = f"CRYST1{box_size[0]:9.3f}{box_size[1]:9.3f}{box_size[2]:9.3f}  90.00  90.00  90.00 P 1           1\n"
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Determine insertion point: after headers and before atoms
    insertion_index = 0
    for i, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM", "MODEL", "END")):  # Place CRYST1 before atom definitions
            insertion_index = i
            break

    # Insert CRYST1 line and rewrite the file
    lines.insert(insertion_index, box_line)
    with open(file_path, "w") as file:
        file.writelines(lines)
