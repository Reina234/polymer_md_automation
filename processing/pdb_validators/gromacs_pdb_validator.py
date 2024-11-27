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
    def validate(self, file_path: str, solvent: Solvent, replace_box: bool = False) -> bool:
        """
        Validate the PDB file, ensuring density alignment and box corrections.
        """
        parser = PDBParser(file_path)
        if not parser.has_box_information():
            logger.warning("Box information missing. Calculating and adding it.")
            box_size = parser.calculate_minimum_box()
            parser.add_or_replace_box_information(box_size)

        density = parser.calculate_density(solvent.molecular_weight)
        expected_density = solvent.density
        deviation = abs((density - expected_density) / expected_density) * 100

        if deviation > DENSITY_TOLERANCE_PERCENTAGE:
            logger.warning(f"Density deviation too high: {deviation:.2f}%")
            if replace_box:
                box_size = parser.adjust_box_for_density(solvent.density, solvent.molecular_weight)
                parser.add_or_replace_box_information(box_size)
                logger.info("Box size adjusted to match density.")
            else:
                logger.error("Density mismatch. Fixing box disabled.")
                return False
        logger.info("Validation successful.")
        return True
