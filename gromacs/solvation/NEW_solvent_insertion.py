import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from config.constants import DEFAULT_BOX_SIZE_NM
from typing import Optional, List
from preprocessing.calculation_utils import calculate_num_particles
from config.constants import LengthUnits
from gromacs.base_gromacs_command import BaseGromacsCommand
import logging
from typing import Tuple

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SolventInsertion(BaseGromacsCommand):
    OUTPUT_NAME = "solvent_box.gro"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def run(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        num_iterations_max: int = 5,
        tolerance: float = 0.05,
        additional_notes: Optional[str] = None,
    ) -> str:
        """
        Iteratively inserts solvent molecules until the desired number of particles
        is achieved or the maximum iterations are reached. Retains previously added
        molecules in subsequent iterations.

        Args:
            solvent_pdb_path (str): Path to the solvent PDB file.
            run_name (str): Name of the run for output files.
            solvent_density (float): Desired solvent density (g/cm³).
            solvent_molecular_weight (float): Solvent molecular weight (g/mol).
            box_size_nm (List[float]): Dimensions of the box (nm).
            num_iterations_max (int): Maximum number of iterations allowed.
            tolerance (float): Acceptable deviation from the desired number of particles.
            additional_notes (Optional[str]): Additional metadata notes.

        Returns:
            str: Path to the final solvent box `.gro` file.
        """
        output_dir = os.path.join(run_name, GROMACS_OUTPUT_SUBDIR)
        os.makedirs(output_dir, exist_ok=True)

        # Calculate target number of particles
        target_num_particles = calculate_num_particles(
            box_size_nm,
            solvent_molecular_weight,
            solvent_density,
            box_units=self.UNITS,
        )
        logger.info(f"Target number of particles: {round(target_num_particles)}")

        # Start with an empty box
        current_gro_path = os.path.join(output_dir, self.OUTPUT_NAME)
        for iteration in range(1, num_iterations_max + 1):
            # Check current particle count
            current_num_particles = (
                self._count_particles(current_gro_path)
                if os.path.exists(current_gro_path)
                else 0
            )
            remaining_particles = round(target_num_particles - current_num_particles)

            logger.info(
                f"Iteration {iteration}: Current particles: {current_num_particles}, Remaining: {remaining_particles}"
            )
            append = os.path.exists(current_gro_path)
            if remaining_particles <= target_num_particles * tolerance:
                logger.info(f"Target achieved with {current_num_particles} particles.")
                if self.metadata_tracker:
                    self._update_metadata(
                        solvent_pdb_path=solvent_pdb_path,
                        run_name=run_name,
                        solvent_density=solvent_density,
                        solvent_molecular_weight=solvent_molecular_weight,
                        box_size_nm=box_size_nm,
                        additional_notes=additional_notes,
                    )
                return current_gro_path

            # Add remaining molecules
            command, _ = self._create_insert_command(
                solvent_pdb_path=solvent_pdb_path,
                box_size_nm=box_size_nm,
                num_particles=remaining_particles,
                output_path=current_gro_path,
                append=append,  # Allow appending to the existing box
            )
            self._execute(command)

        # If the loop ends without achieving the target
        raise RuntimeError(
            f"Failed to reach the target number of particles after {num_iterations_max} iterations."
        )

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
            "action(s)": f"Inserted solvent molecules from {solvent_pdb_path} to {run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.OUTPUT_NAME}",
            "additional_notes": additional_notes,
        }
