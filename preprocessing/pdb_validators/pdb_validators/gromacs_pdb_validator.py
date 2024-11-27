import os
import logging
from typing import List, Optional
from preprocessing.file_commenter import FileCommenter
from preprocessing.utils import check_file_exists, check_file_type
from abc import ABC, abstractmethod
from preprocessing.pdb_validators.pdb_validators.base_pdb_validator import (
    BasePDBValidator,
)
from preprocessing.pdb_validators.pdb_utils import (
    calculate_minimum_box_size,
    calculate_volume_for_desired_density,
    scale_box_to_desired_volume,
    calculate_density,
)
from data_models.solvent import Solvent
from preprocessing.pdb_validators.pdb_parser import PDBParser
from config.constants import DENSITY_TOLERANCE_PERCENTAGE
from config.constants import LengthUnits

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class GROMACSPDBValidator(BasePDBValidator):
    BOX_GENERATED_COMMENT = "Box dimensions generated based on atom coordinates."
    DENSITY_VALIDATION_PASSED_COMMENT = "Density validation passed."
    DENSITY_VALIDATION_FAILED_COMMENT = "Density validation failed."

    def __init__(self, auto_fix_dimensions: bool = False):
        super().__init__()
        self.auto_fix_dimensions = auto_fix_dimensions

    def validate(
        self,
        input_file_path: str,
        solvent: Solvent,
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
            self._validate_gromacs_specifics(input_file_path, solvent, output_file_path)

        except ValueError as e:
            logger.error(f"[!] Validation halted: {e}")
            raise

    def _validate_gromacs_specifics(
        self,
        input_file_path: str,
        solvent: Solvent,
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
            logger.info(f"No box_dimensions found in the PDB file.")
            box_dimensions, content = self._generate_box_dimensions(solvent, content)
            self.pdb_parser.save(output_file_path, content)

        else:
            logger.info(f"Box dimensions found: {box_dimensions}")
            density_check = self._density_check(solvent, box_dimensions)
            if not density_check and self.auto_fix_dimensions:
                content = self.pdb_parser.add_comment(
                    content, self.DENSITY_VALIDATION_FAILED_COMMENT
                )
                logger.info(f"Attempting to fix box dimensions...")
                box_dimensions, content = self._generate_box_dimensions(
                    solvent, content
                )
                self.pdb_parser.save(output_file_path, content)

            elif not density_check:
                content = self.pdb_parser.add_comment(
                    content, self.DENSITY_VALIDATION_FAILED_COMMENT
                )
                self.pdb_parser.save(output_file_path, content)
                raise ValueError(
                    f"Validation failed: Density check failed for {input_file_path}."
                )
            elif density_check:
                content = self.pdb_parser.add_comment(
                    content, self.DENSITY_VALIDATION_PASSED_COMMENT
                )
                self.pdb_parser.save(output_file_path, content)

    def _generate_box_dimensions(
        self, solvent: Solvent, content: List[str]
    ) -> List[float]:
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
        min_box_size = calculate_minimum_box_size(atom_coordinates)
        volume_needed = calculate_volume_for_desired_density(
            solvent.molecular_weight, solvent.density
        )
        box_dimensions = scale_box_to_desired_volume(min_box_size, volume_needed)
        content = self.pdb_parser.add_or_replace_box_dimensions(content, box_dimensions)
        content = self.pdb_parser.add_comment(content, self.BOX_GENERATED_COMMENT)
        return box_dimensions, content

    def _density_check(self, solvent: Solvent, box_dimensions: List[float]) -> None:
        """
        Perform density check based on box dimensions.

        Args:
            box_dimensions (List[float]): Box dimensions [x, y, z] in nm.

        Raises:
            ValueError: If the calculated density does not match the desired density.
        """
        logger.info(f"Performing density check...")
        density = calculate_density(
            solvent.molecular_weight, box_dimensions, length_units=LengthUnits.ANGSTROM
        )
        if (
            density - solvent.density
        ) ** 2 > DENSITY_TOLERANCE_PERCENTAGE * solvent.density / 100:
            logger.error(
                f"Calculated density ({density:.6f} kg/m³) does not match desired density ({solvent.density} g/cm³)."
            )
            return False
        logger.info(f"Density check passed: {density:.6f} kg/m³")
        return True
