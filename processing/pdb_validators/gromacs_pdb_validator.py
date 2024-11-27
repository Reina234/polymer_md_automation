from typing import Optional
from processing.pdb_validators.base_pdb_validator import BasePDBValidator
from processing.pdb_validators.pdb_utils import add_box_information, calculate_box_size
from data_models.solvent import Solvent
from processing.metadata_tracker import MetadataTracker
import logging
from processing.file_operations import move_file
from processing.parsers.EDIT_pdb_parser import PDBParser
from processing.file_commenter import FileCommenter
from config.constants import DENSITY_TOLERANCE_PERCENTAGE
logger = logging.getLogger(__name__)

#NOTE: consider adding renaming functionality 


class GROMACSPDBValidator(BasePDBValidator):
    """
    Validator for GROMACS-compatible PDB files, including density checks and fixing.
    """

    def __init__(self, solvent: Solvent, metadata_tracker: Optional[object] = None):
        super().__init__(metadata_tracker)
        self.solvent = solvent
        self.commenter = FileCommenter()

    def validate(
        self,
        file_path: str,
        output_location: Optional[str] = None,
        move_original: Optional[str] = None,
        replace_box: bool = False,
    ) -> bool:
        """
        Validate the PDB file for GROMACS compatibility.

        Args:
            file_path (str): Path to the PDB file.
            output_location (Optional[str]): Directory to save the output file (default is same location).
            move_original (Optional[str]): Directory to move the original file.
            replace_box (bool): Whether to replace the box information if density is out of tolerance.

        Returns:
            bool: True if the file is valid or fixed, False otherwise.
        """
        try:
            parser = PDBParser(file_path)
        except (FileNotFoundError, ValueError) as e:
            logger.error(f"[!] {e}")
            return False

        try:
            if not parser.has_box_information():
                logger.warning(f"[!] Box information missing in {file_path}.")
                box_size = calculate_box_size(self.solvent.molecular_weight, self.solvent.density)
                parser.add_or_replace_box_information(box_size)
                parser.add_comment("Box information added.")
            else:
                # Check for zero dimensions
                x, y, z = parser.box_dimensions
                if x == 0 or y == 0 or z == 0:
                    raise ValueError(f"Box dimensions cannot be zero. Found: x={x}, y={y}, z={z}")

            density = parser.calculate_density(self.solvent.molecular_weight)
            expected_density = self.solvent.density
            deviation = abs((density - expected_density) / expected_density) * 100

            if deviation > DENSITY_TOLERANCE_PERCENTAGE:
                logger.warning(f"[!] Density tolerance exceeded for {file_path}. Deviation: {deviation:.2f}%")
                if replace_box:
                    box_size = calculate_box_size(self.solvent.molecular_weight, self.solvent.density)
                    parser.add_or_replace_box_information(box_size)
                    parser.add_comment("Box information replaced due to density deviation.")
                    density = parser.calculate_density(self.solvent.molecular_weight)
                    deviation = abs((density - expected_density) / expected_density) * 100
                    logger.info(f"[+] New density: {density:.6f} g/cm³, Deviation: {deviation:.2f}%")
                else:
                    logger.warning(f"[!] Density tolerance exceeded for {file_path}. File not fixed.")
                    parser.add_comment(f"Density tolerance exceeded ({deviation:.2f}%). Box not replaced.")
                    return False
            else:
                logger.info(f"[+] Density check passed for {file_path}. Deviation: {deviation:.2f}%")
                parser.add_comment(f"Density check passed. Deviation: {deviation:.2f}%")

            if move_original:
                move_file(file_path, move_original)
                parser.add_comment(f"Original file moved to: {move_original}")

            output_file = parser.save(output_location)
            parser.add_comment(f"File saved to: {output_location}" if output_location else "File overwritten.")
            logger.info(f"[+] GROMACS validation completed for {output_file}.")
            return output_file

        except ValueError as e:
            logger.error(f"[!] Validation failed: {e}")
            return None 