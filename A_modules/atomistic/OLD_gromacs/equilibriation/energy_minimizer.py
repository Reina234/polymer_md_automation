import os
import shutil
from typing import Optional, Tuple, List, Dict
from config.paths import (
    TEMPORARY_OUTPUT_DIR,
    MDP_FULL_PATHS,
    TemplatedMdps,
    GROMACS_OUTPUT_SUBDIR,
    EQUILIBRIUM_SUBDIR,
    BASE_OUTPUT_DIR,
)
from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.shared.command_line_operation import CommandLineOperation
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class EnergyMinimizer(CommandLineOperation):
    output_name = "em.gro"
    output_tpr_name = "em.tpr"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def run(
        self,
        input_gro: str,
        input_top: str,
        run_name: str,
    ) -> str:
        """
        Run energy minimization for the given input structure and topology.

        Args:
            input_gro (str): Path to the input .gro file (e.g., solvated structure).
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run (used to organize outputs).

        Returns:
            str: Path to the minimized structure (`em.gro`) in the final output directory.
        """
        # Paths for temporary and final outputs
        minim_mdp_path = MDP_FULL_PATHS[TemplatedMdps.MINIM.value]
        temp_output_dir = TEMPORARY_OUTPUT_DIR
        final_output_dir = os.path.join(
            run_name, GROMACS_OUTPUT_SUBDIR, EQUILIBRIUM_SUBDIR
        )

        # Create commands for GROMACS
        grompp_command, output_tpr_path, temp_output_dir = self._create_grompp_command(
            input_gro=input_gro,
            input_top=input_top,
            minim_mdp_path=minim_mdp_path,
        )
        self._execute(grompp_command)

        mdrun_command, em_gro_temp_path = self._create_mdrun_command(
            output_tpr_path=output_tpr_path, output_dir=temp_output_dir
        )
        self._execute(mdrun_command)

        # Check if the output exists in the temporary directory
        if not os.path.exists(em_gro_temp_path):
            raise FileNotFoundError(
                f"Energy minimization failed. {em_gro_temp_path} not found."
            )

        # Move `em.gro` to the final output directory
        os.makedirs(final_output_dir, exist_ok=True)
        em_gro_final_path = os.path.join(final_output_dir, self.OUTPUT_NAME)
        shutil.move(em_gro_temp_path, em_gro_final_path)
        logger.info(f"Moved minimized structure to {em_gro_final_path}")

        # Update metadata
        if self.metadata_tracker:
            self._update_metadata(
                input_gro=input_gro,
                input_top=input_top,
                run_name=run_name,
                minim_mdp_path=minim_mdp_path,
            )

        return em_gro_final_path

    def metadata(
        self,
        input_gro: str,
        input_top: str,
        run_name: str,
        minim_mdp_path: Optional[str] = None,
    ) -> Dict[str, any]:
        """
        Generate metadata for the energy minimization process.

        Args:
            input_gro (str): Path to the input .gro file.
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run.
            minim_mdp_path (Optional[str]): Path to the energy minimization .mdp file.

        Returns:
            dict: Metadata describing the energy minimization process.
        """
        return {
            "program_used": "GROMACS",
            "commands_executed": [
                {
                    "type": "grompp",
                    "description": "Preparation of input files for minimization",
                    "mdp_file": minim_mdp_path,
                    "input_structure": input_gro,
                    "input_topology": input_top,
                },
                {
                    "type": "mdrun",
                    "description": "Energy minimization",
                    "input_tpr": os.path.join(
                        TEMPORARY_OUTPUT_DIR, self.output_tpr_name
                    ),
                    "output_structure": self.OUTPUT_NAME,
                },
            ],
            "inputs": {
                "structure": input_gro,
                "topology": input_top,
                "mdp_file": minim_mdp_path,
            },
            "outputs": {
                "final_structure": os.path.join(
                    BASE_OUTPUT_DIR,
                    run_name,
                    GROMACS_OUTPUT_SUBDIR,
                    EQUILIBRIUM_SUBDIR,
                    self.OUTPUT_NAME,
                ),
                "temporary_files": TEMPORARY_OUTPUT_DIR,
            },
        }

    def _create_grompp_command(
        self,
        input_gro: str,
        input_top: str,
        minim_mdp_path: str,
    ) -> Tuple[List[str], str, str]:
        """
        Create the grompp command for preparing the minimization run.

        Args:
            input_gro (str): Path to the input .gro file.
            input_top: Path to the input topology file.
            minim_mdp_path: Path to the energy minimization .mdp file.

        Returns:
            Tuple[List[str], str, str]: grompp command, output .tpr path, and output directory.
        """
        temp_output_dir = TEMPORARY_OUTPUT_DIR
        os.makedirs(temp_output_dir, exist_ok=True)

        output_tpr_path = os.path.join(temp_output_dir, self.output_tpr_name)

        grompp_command = [
            "gmx",
            "grompp",
            "-f",
            minim_mdp_path,
            "-c",
            input_gro,
            "-p",
            input_top,
            "-o",
            output_tpr_path,
            "-maxwarn",
            "1",
        ]
        return grompp_command, output_tpr_path, temp_output_dir

    def _create_mdrun_command(
        self, output_tpr_path: str, output_dir: str
    ) -> Tuple[List[str], str]:
        """
        Create the mdrun command for performing the minimization.

        Args:
            output_tpr_path (str): Path to the prepared .tpr file.
            output_dir (str): Directory to save the temporary output.

        Returns:
            Tuple[List[str], str]: mdrun command and output .gro path.
        """
        em_gro_path = os.path.join(output_dir, self.OUTPUT_NAME)
        mdrun_command = [
            "gmx",
            "mdrun",
            "-deffnm",
            os.path.join(output_dir, "em"),
        ]
        return mdrun_command, em_gro_path
