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
from gromacs.base_gromacs_command import BaseGromacsCommand
import logging
from preprocessing.template_utils import retrieve_mdps

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
from gromacs.base_gromacs_command import BaseGromacsCommand
import logging
from preprocessing.template_utils import retrieve_mdps

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class TemperatureEquilibration(BaseGromacsCommand):
    OUTPUT_NAME = "nvt.gro"
    OUTPUT_TPR_NAME = "nvt.tpr"

    def __init__(self, metadata_tracker=None):
        super().__init__(metadata_tracker)

    def run(
        self,
        em_gro: str,
        input_top: str,
        run_name: str,
        simulation_temp_k: float,
        output_base_dir: str = BASE_OUTPUT_DIR,
        rdd: Optional[float] = None,
        ntmpi: Optional[int] = None,
    ) -> str:
        """
        Run NVT equilibration for the given input structure and topology.

        Args:
            em_gro (str): Path to the input structure file (output of energy minimization).
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run (used to organize outputs).
            simulation_temp_k (float): Target temperature for the simulation (in K).
            output_base_dir (str): Base directory for output files.
            rdd (Optional[float]): Minimum allowed cell size for domain decomposition.
            ntmpi (Optional[int]): Number of MPI threads to use.

        Returns:
            str: Path to the equilibrated structure (`nvt.gro`) in the final output directory.
        """
        # Paths for temporary and final outputs
        nvt_mdp_path = retrieve_mdps(TemplatedMdps.NVT, simulation_temp_k)
        temp_output_dir = TEMPORARY_OUTPUT_DIR
        final_output_dir = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, EQUILIBRIUM_SUBDIR
        )

        # Create commands for GROMACS
        grompp_command, output_tpr_path, temp_output_dir = self._create_grompp_command(
            input_gro=em_gro,
            input_top=input_top,
            nvt_mdp_path=nvt_mdp_path,
        )
        self._execute(grompp_command)

        mdrun_command, nvt_gro_temp_path = self._create_mdrun_command(
            output_tpr_path=output_tpr_path,
            output_dir=temp_output_dir,
            rdd=rdd,
            ntmpi=ntmpi,
        )
        self._execute(mdrun_command)

        # Check if the output exists in the temporary directory
        if not os.path.exists(nvt_gro_temp_path):
            raise FileNotFoundError(
                f"NVT equilibration failed. {nvt_gro_temp_path} not found."
            )

        # Move `nvt.gro` to the final output directory
        os.makedirs(final_output_dir, exist_ok=True)
        nvt_gro_final_path = os.path.join(final_output_dir, self.OUTPUT_NAME)
        shutil.move(nvt_gro_temp_path, nvt_gro_final_path)
        logger.info(f"Moved NVT equilibrated structure to {nvt_gro_final_path}")

        # Update metadata
        if self.metadata_tracker:
            self._update_metadata(
                em_gro=em_gro,
                input_top=input_top,
                run_name=run_name,
                nvt_mdp_path=nvt_mdp_path,
            )

        return nvt_gro_final_path

    def _create_grompp_command(
        self,
        input_gro: str,
        input_top: str,
        nvt_mdp_path: str,
    ) -> Tuple[List[str], str, str]:
        """
        Create the grompp command for preparing the NVT equilibration run.

        Args:
            input_gro (str): Path to the input structure file.
            input_top: Path to the input topology file.
            nvt_mdp_path: Path to the NVT equilibration .mdp file.

        Returns:
            Tuple[List[str], str, str]: grompp command, output .tpr path, and output directory.
        """
        temp_output_dir = TEMPORARY_OUTPUT_DIR
        os.makedirs(temp_output_dir, exist_ok=True)

        output_tpr_path = os.path.join(temp_output_dir, self.OUTPUT_TPR_NAME)

        grompp_command = [
            "gmx",
            "grompp",
            "-f",
            nvt_mdp_path,
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
        self,
        output_tpr_path: str,
        output_dir: str,
        rdd: Optional[float],
        ntmpi: Optional[int],
    ) -> Tuple[List[str], str]:
        """
        Create the mdrun command for performing the NVT equilibration.

        Args:
            output_tpr_path (str): Path to the prepared .tpr file.
            output_dir (str): Directory to save the temporary output.
            rdd (Optional[float]): Minimum allowed cell size for domain decomposition.
            ntmpi (Optional[int]): Number of MPI threads.

        Returns:
            Tuple[List[str], str]: mdrun command and output .gro path.
        """
        nvt_gro_path = os.path.join(output_dir, self.OUTPUT_NAME)
        mdrun_command = [
            "gmx",
            "mdrun",
            "-deffnm",
            os.path.join(output_dir, "nvt"),
        ]

        if rdd:
            mdrun_command.extend(["-rdd", str(rdd)])
        if ntmpi:
            mdrun_command.extend(["-ntmpi", str(ntmpi)])

        return mdrun_command, nvt_gro_path

    def metadata(
        self,
        em_gro: str,
        input_top: str,
        run_name: str,
        nvt_mdp_path: Optional[str] = None,
    ) -> Dict[str, any]:
        """
        Generate metadata for the NVT equilibration process.

        Args:
            em_gro (str): Path to the input structure file (output of energy minimization).
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run.
            nvt_mdp_path (Optional[str]): Path to the NVT .mdp file.

        Returns:
            dict: Metadata describing the NVT equilibration process.
        """
        return {
            "program_used": "GROMACS",
            "commands_executed": [
                {
                    "type": "grompp",
                    "description": "Preparation of input files for NVT equilibration",
                    "mdp_file": nvt_mdp_path,
                    "input_structure": em_gro,
                    "input_topology": input_top,
                },
                {
                    "type": "mdrun",
                    "description": "NVT equilibration",
                    "input_tpr": os.path.join(
                        TEMPORARY_OUTPUT_DIR, self.OUTPUT_TPR_NAME
                    ),
                    "output_structure": self.OUTPUT_NAME,
                },
            ],
            "inputs": {
                "structure": em_gro,
                "topology": input_top,
                "mdp_file": nvt_mdp_path,
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
