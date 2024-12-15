import os

# from A_config.paths import GROMACS_OUTPUT_SUBDIR
from A_config.constants import DEFAULT_BOX_SIZE_NM
from typing import Optional, List
from A_modules.shared.utils.calculation_utils import calculate_num_particles
from A_modules.shared.utils.utils import directory_exists_check_wrapper
from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.shared.command_line_operation import CommandLineOperation
import logging
from typing import Tuple
import subprocess

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SolventInsertion(CommandLineOperation):
    output_name = "solvent_box.gro"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    @directory_exists_check_wrapper(dir_arg_index=2)
    def run(
        self,
        input_pdb_path: str,
        output_dir: str,
        desired_density: float,
        molecular_weight: float,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        num_iterations_max: int = 5,
        tolerance: float = 0.05,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ) -> str:
        self._validate_parameters(
            desired_density, molecular_weight, num_iterations_max, tolerance
        )

        gro_path = os.path.join(output_dir, self.output_name)

        target_num_particles = self._calculate_target_particles(
            box_size_nm, molecular_weight, desired_density, verbose
        )

        for iteration in range(1, num_iterations_max + 1):
            if self._is_target_achieved(
                gro_path, target_num_particles, tolerance, iteration, verbose
            ):
                # Update metadata only on success
                if self.metadata_tracker:
                    self._update_metadata(
                        input_pdb_path=input_pdb_path,
                        output_dir=output_dir,
                        desired_density=desired_density,
                        molecular_weight=molecular_weight,
                        iteration=iteration,
                        box_size_nm=box_size_nm,
                        additional_notes=additional_notes,
                    )
                return gro_path

            remaining_particles = self._calculate_remaining_particles(
                gro_path, target_num_particles
            )
            self._add_particles(
                input_pdb_path, gro_path, box_size_nm, remaining_particles, verbose
            )

        # Raise an error if target is not met after max iterations
        raise RuntimeError(
            f"Failed to reach the target number of particles after {num_iterations_max} iterations."
        )

    def metadata(
        self,
        input_pdb_path: str,
        output_dir: str,
        desired_density: float,
        molecular_weight: float,
        iteration: int,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        additional_notes: Optional[str] = None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS insert-molecules",
            "details": f"Created solvent box of size {box_size_nm} (nm) with target density {desired_density} g/cm³ after {iteration} iterations",
            "action(s)": f"Inserted solvent molecules from {input_pdb_path} to {output_dir}/{self.output_name}",
            "additional_notes": additional_notes,
        }

    def _calculate_target_particles(
        self, box_size_nm, molecular_weight, desired_density, verbose
    ):
        """
        Calculates the target number of particles based on box size, molecular weight, and desired density.
        """
        target_num_particles = calculate_num_particles(
            box_size_nm,
            molecular_weight=molecular_weight,
            density=desired_density,
            box_units=self.default_units,
        )
        if verbose:
            logger.info(f"Target number of particles: {round(target_num_particles)}")
        return target_num_particles

    def _is_target_achieved(
        self, gro_path, target_num_particles, tolerance, iteration, verbose
    ):
        """
        Checks if the target number of particles has been achieved.
        """
        current_num_particles = (
            self._count_particles(gro_path) if os.path.exists(gro_path) else 0
        )
        remaining_particles = target_num_particles - current_num_particles

        if verbose:
            logger.info(
                f"Iteration {iteration}: Current particles: {current_num_particles}, Remaining: {remaining_particles}"
            )

        return remaining_particles <= target_num_particles * tolerance

    def _calculate_remaining_particles(self, gro_path, target_num_particles):
        """
        Calculates the number of remaining particles needed to reach the target.
        """
        current_num_particles = (
            self._count_particles(gro_path) if os.path.exists(gro_path) else 0
        )
        return target_num_particles - current_num_particles

    def _add_particles(
        self, input_pdb_path, gro_path, box_size_nm, remaining_particles, verbose
    ):
        """
        Adds remaining particles to the system by executing the insert command.
        """
        command, _ = self._create_insert_command(
            solvent_pdb_path=input_pdb_path,
            box_size_nm=box_size_nm,
            num_particles=remaining_particles,
            output_path=gro_path,
            append=True,  # Allow appending to the existing box
        )
        self._execute(command, verbose=verbose)

    def _validate_parameters(
        self,
        desired_density: float,
        molecular_weight: float,
        num_iterations_max: int,
        tolerance: float,
    ):
        if desired_density <= 0:
            raise ValueError("Desired density must be positive.")
        if molecular_weight <= 0:
            raise ValueError("Molecular weight must be positive.")
        if num_iterations_max <= 0:
            raise ValueError("Number of iterations must be positive.")
        if not (0 < tolerance < 1):
            raise ValueError("Tolerance must be between 0 and 1.")

    def _create_insert_command(
        self,
        solvent_pdb_path: str,
        box_size_nm: List[float],
        num_particles: int,
        output_path: str,
        append: bool = False,
    ) -> Tuple[List[str], str]:
        command = [
            "gmx",
            "insert-molecules",
            "-ci",
            solvent_pdb_path,
            "-nmol",
            str(num_particles),
            "-box",
            str(box_size_nm[0]),
            str(box_size_nm[1]),
            str(box_size_nm[2]),
            "-o",
            output_path,
        ]
        if append:
            command.append("-f")
            command.append(output_path)  # Append to the existing file
        return command, output_path

    def _count_particles(self, gro_file: str) -> int:
        """
        Count the number of molecules (particles) in a .gro file.

        Args:
            gro_file (str): Path to the `.gro` file.

        Returns:
            int: Number of molecules in the file.
        """
        if not os.path.exists(gro_file):
            return 0

        with open(gro_file, "r") as file:
            lines = file.readlines()

        # Second line of .gro file contains the total number of ATOMS
        try:
            num_atoms = int(lines[1].strip())
        except (IndexError, ValueError) as e:
            raise ValueError(f"Error reading atom count from {gro_file}: {e}")

        # Count unique residue/molecule IDs (first 5 characters of each line with atomic data)
        residue_ids = set()
        for line in lines[2 : 2 + num_atoms]:  # Only iterate over atomic data lines
            residue_id = line[:5].strip()
            residue_ids.add(residue_id)

        return len(residue_ids)

    def metadata(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        additional_notes: Optional[str] = None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS insert-molecules",
            "details": f"Created solvent box of size {box_size_nm} (nm) with target density {solvent_density} g/cm³",
            "action(s)": f"Inserted solvent molecules from {solvent_pdb_path} to {run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.output_name}",
            "additional_notes": additional_notes,
        }
