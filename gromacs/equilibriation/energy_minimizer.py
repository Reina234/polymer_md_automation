import os
from typing import Optional, Tuple, List
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
from preprocessing.template_utils import create_mdps, retrieve_mdps

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class EnergyMinimizer(BaseGromacsCommand):
    OUTPUT_NAME = "em.gro"
    OUTPUT_TPR_NAME = "em.tpr"

    def __init__(self, metadata_tracker=None):
        super().__init__(metadata_tracker)

    def run(
        self,
        input_gro: str,
        input_top: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
    ) -> str:
        """
        Run energy minimization for the given input structure and topology.

        Args:
            input_gro (str): Path to the input .gro file (e.g., solvated structure).
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run (used to organize outputs).
            minim_mdp_path (Optional[str]): Path to the energy minimization .mdp file.
                                             If not provided, defaults to the one in MDP_FULL_PATHS.

        Returns:
            str: Path to the minimized structure (`em.gro`).
        """
        minim_mdp_path = MDP_FULL_PATHS[TemplatedMdps.MINIM.value]
        output_dir = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, EQUILIBRIUM_SUBDIR
        )
        grompp_command, output_tpr_path, output_dir = self._create_grompp_command(
            input_gro=input_gro,
            input_top=input_top,
            run_name=run_name,
            minim_mdp_path=minim_mdp_path,
        )
        self._execute(grompp_command)

        mdrun_command, em_gro_path = self._create_mdrun_command(
            output_tpr_path=output_tpr_path, output_dir=output_dir
        )
        self._execute(mdrun_command)

        if self.metadata_tracker:
            self._update_metadata(
                input_gro=input_gro,
                input_top=input_top,
                run_name=run_name,
                minim_mdp_path=minim_mdp_path,
            )

        if not os.path.exists(em_gro_path):
            raise FileNotFoundError(
                f"Energy minimization failed. {em_gro_path} not found."
            )
        return em_gro_path

    def _create_grompp_command(
        self,
        input_gro: str,
        input_top: str,
        run_name: str,
        minim_mdp_path: Optional[str] = None,
    ) -> Tuple[List[str], str, str]:
        """
        Create the grompp command for preparing the minimization run.

        Args:
            input_gro (str): Path to the input .gro file.
            input_top (str): Path to the input topology file.
            run_name (str): Name of the run.
            minim_mdp_path (Optional[str]): Path to the energy minimization .mdp file.

        Returns:
            Tuple[List[str], str, str]: grompp command, output .tpr path, and output directory.
        """
        output_dir = TEMPORARY_OUTPUT_DIR
        os.makedirs(output_dir, exist_ok=True)

        if not minim_mdp_path:
            minim_mdp_path = MDP_FULL_PATHS[TemplatedMdps.MINIM.value]

        output_tpr_path = os.path.join(output_dir, self.OUTPUT_TPR_NAME)

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
        return grompp_command, output_tpr_path, output_dir

    def _create_mdrun_command(
        self, output_tpr_path: str, output_dir: str
    ) -> Tuple[List[str], str]:
        """
        Create the mdrun command for performing the minimization.

        Args:
            output_tpr_path (str): Path to the prepared .tpr file.
            output_dir (str): Directory to save the output.

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

    def metadata(
        self,
        input_gro: str,
        input_top: str,
        run_name: str,
        minim_mdp_path: Optional[str] = None,
    ) -> dict:
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
            "program(s) used": "GROMACS grompp and mdrun",
            "details": f"Energy minimization performed using {minim_mdp_path}",
            "inputs": {
                "structure": input_gro,
                "topology": input_top,
            },
            "output": f"{TEMPORARY_OUTPUT_DIR}/{run_name}/{self.OUTPUT_NAME}",
        }
