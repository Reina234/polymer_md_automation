from typing import Dict, List, Optional
import logging
import os

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class GROParser:
    """
    A stateless class to parse GRO files and extract relevant information.
    """

    @staticmethod
    def read_file(input_file_path: str) -> List[str]:
        """
        Read the GRO file and return its content as a list of lines.

        Args:
            input_file_path (str): Path to the GRO file.

        Returns:
            List[str]: Lines from the GRO file.
        """
        if not os.path.exists(input_file_path):
            raise FileNotFoundError(f"File not found: {input_file_path}")
        with open(input_file_path, "r") as file:
            content = file.readlines()
        return content

    @staticmethod
    def extract_molecule_counts(content: List[str]) -> Dict[str, int]:
        """
        Extract the number of molecules for each residue type from the GRO file.

        Args:
            content (List[str]): Lines from the GRO file.

        Returns:
            Dict[str, int]: A dictionary mapping residue names to their molecule counts.
        """
        molecule_counts = {}
        current_residue = None

        for line in content[2:-1]:  # Skip title, atom count, and box dimensions
            if len(line) < 20:
                continue  # Skip malformed lines
            residue_name = line[:5].strip()
            if residue_name != current_residue:
                current_residue = residue_name
                molecule_counts[residue_name] = molecule_counts.get(residue_name, 0) + 1

        return molecule_counts

    @staticmethod
    def extract_box_dimensions(content: List[str]) -> Optional[List[float]]:
        """
        Extract box dimensions from the last line of the GRO file.

        Args:
            content (List[str]): Lines from the GRO file.

        Returns:
            Optional[List[float]]: Box dimensions [x, y, z] in nm, or None if not found.
        """
        try:
            box_dimensions = list(map(float, content[-1].split()))
            if len(box_dimensions) >= 3:
                return box_dimensions[:3]
            else:
                logger.warning("[!] Box dimensions are incomplete in the GRO file.")
                return None
        except ValueError:
            logger.warning("[!] Invalid box dimensions in the GRO file.")
            return None

    @staticmethod
    def validate_molecule_counts(
        molecule_counts: Dict[str, int], expected_molecule_counts: Dict[str, int]
    ):
        """
        Validate the extracted molecule counts against the expected counts.

        Args:
            molecule_counts (Dict[str, int]): Extracted molecule counts from the GRO file.
            expected_molecule_counts (Dict[str, int]): Expected molecule counts for validation.

        Raises:
            ValueError: If the molecule counts do not match expectations.
        """
        for residue, count in molecule_counts.items():
            if residue in expected_molecule_counts:
                if count != expected_molecule_counts[residue]:
                    raise ValueError(
                        f"[!] Mismatch for residue {residue}: Expected {expected_molecule_counts[residue]} molecules, "
                        f"but found {count}."
                    )
            else:
                logger.warning(f"[!] Unexpected residue {residue} in GRO file.")

        logger.info("[+] Molecule counts validated successfully.")
