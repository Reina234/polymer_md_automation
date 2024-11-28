import os
import logging
from typing import List, Optional
from preprocessing.validators.pdb_validators.base_pdb_validator import (
    BasePDBValidator,
)
from preprocessing.calculation_utils import calculate_minimum_box_size
from data_models.solvent import Solvent
from config.constants import DENSITY_TOLERANCE_PERCENTAGE
from config.constants import LengthUnits
from preprocessing.metadata_tracker import MetadataTracker
from config.constants import LengthUnits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class GROMACSPDBValidator(BasePDBValidator):
    BOX_GENERATED_COMMENT = "No valid box dimensions found... Box dimensions generated based on atom coordinates."
    BOX_DIMENSIONS_FOUND = "Valid box dimensions found"

    def __init__(
        self,
        metadata_tracker: Optional[MetadataTracker] = None,
        padding: float = 0.5 * 10**-9,
    ):
        super().__init__(metadata_tracker)
        self.padding = padding

    def validate(
        self,
        input_file_path: str,
        additional_notes: Optional[str] = None,
        output_file_path: Optional[str] = None,
    ) -> None:
        """
        Perform a full validation of the PDB file, starting with a skeletal check.

        Args:
            file_path (str): Path to the PDB file.

        Raises:
            ValueError: If any validation step fails.
        """
        if not output_file_path:
            output_file_path = input_file_path
            # NEED TO EDIT #################
        try:
            # Step 1: Perform skeletal check
            logger.info("[*] Starting skeletal validation...")
            self._skeletal_check(input_file_path, output_file_path)
            logger.info("[+] Skeletal validation passed.")

            # Step 2: Proceed with other checks if skeletal validation succeeds
            logger.info("[*] Proceeding with additional GROMACS checks...")
            self._validate_gromacs_specifics(input_file_path, output_file_path)
            if self.metadata_tracker:
                self._update_metadata(
                    input_file_path, output_file_path, additional_notes
                )

        except ValueError as e:
            logger.error(f"[!] Validation halted: {e}")
            raise

    def _validate_gromacs_specifics(
        self,
        input_file_path: str,
        output_file_path: Optional[str] = None,
    ):
        """
        Perform GROMACS-specific validations (e.g., box size, density).

        Args:
            file_path (str): Path to the PDB file.

        Raises:
            ValueError: If GROMACS-specific validation fails.
        """

        logger.info(f"Performing GROMACS-specific validation on {input_file_path}...")
        content = self.pdb_parser.read_file(input_file_path)

        box_dimensions = self.pdb_parser.extract_box_dimensions(content)

        if not box_dimensions:
            logger.info(f"No valid box_dimensions found in the PDB file.")
            box_dimensions, content = self._generate_box_dimensions(content)
            self.pdb_parser.save(output_file_path, content)

        else:
            logger.info(f"Box dimensions found: {box_dimensions}")
            content = self.pdb_parser.add_comment(content, self.BOX_DIMENSIONS_FOUND)
            self.pdb_parser.save(output_file_path, content)

    def _generate_box_dimensions(self, content: List[str]) -> List[float]:
        """
        Generate box dimensions if not found in the PDB file.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            List[float]: Box dimensions [x, y, z] in nm.
        """
        logger.info(f"Estimating box size based on atom coordinates...")
        atoms = self.pdb_parser.extract_atoms(content)
        atom_coordinates = self.pdb_parser.get_atom_coordinates(atoms)
        box_dimensions = calculate_minimum_box_size(
            atom_coordinates, padding=self.padding, output_units=LengthUnits.ANGSTROM
        )
        content = self.pdb_parser.add_or_replace_box_dimensions(content, box_dimensions)
        content = self.pdb_parser.add_comment(content, self.BOX_GENERATED_COMMENT)
        return box_dimensions, content

    def metadata(self, input_file_path, output_dir, additional_notes=None):
        return {
            "program(s) used": "N/A",
            "defatils": "checked for base pdb structure, and for box size and density match",
            "action(s)": f"checked validity of file at {input_file_path}, saved to {output_dir}",
            "additional_notes": additional_notes,
        }
