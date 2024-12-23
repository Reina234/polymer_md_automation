from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from typing import Optional, Dict, List
from abc import ABC, abstractmethod


from typing import Dict, List, Optional
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.shared.utils.file_utils import batch_copy_file, check_directory_exists
import logging
import os

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class BaseWorkflowStep:
    def __init__(
        self,
        grompp: Grompp,
        mdrun: MDrun,
        default_params: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize the workflow step.

        :param grompp: Grompp command instance.
        :param mdrun: MDrun command instance.
        :param mdp_cache: MDP cache for managing and reusing MDP files.
        :param default_params: Default parameters for this workflow step.
        """
        self.grompp = grompp
        self.mdrun = mdrun  # Dependency Injection
        self.default_params = default_params or {}

    def merge_params(self, varying_params: Dict[str, str]) -> Dict[str, str]:
        """
        Merge varying parameters with default parameters.
        Priority is given to varying parameters.

        :param varying_params: Parameters to override the defaults.
        :return: Merged parameters.
        """
        return {**self.default_params, **varying_params}

    def run(
        self,
        mdp_template_path: str,
        input_gro_path: str,
        input_topol_path: str,
        temp_output_dir: str,
        main_output_dir: str,
        keep_files: List[str],
        varying_params: Dict[str, str],
        mdp_cache: MDPCache,
        verbose: bool = False,
    ) -> None:
        """
        Run the workflow step.

        :param mdp_template_path: Path to the MDP template file.
        :param input_gro_path: Path to the input GRO file.
        :param input_topol_path: Path to the topology file.
        :param temp_output_dir: Directory for temporary output files.
        :param main_output_dir: Directory for final output files.
        :param keep_files: List of file extensions to keep in the main output directory.
        :param varying_params: Parameters specific to this workflow step.
        :param verbose: Enable verbose logging for GROMACS commands.
        """
        # Merge parameters
        params = self.merge_params(varying_params)

        # Ensure directories exist
        check_directory_exists(temp_output_dir)
        check_directory_exists(main_output_dir)

        mdp_file = mdp_cache.get_or_create_mdp(
            template_path=mdp_template_path, params=params
        )

        # Generate output file paths
        # tpr_file = os.path.join(temp_output_dir, "step.tpr")
        output_prefix = os.path.join(temp_output_dir, "step_output")

        # Run GROMPP
        grompp_output = self.grompp.run(
            mdp_file_path=mdp_file,
            input_gro_path=input_gro_path,
            input_topol_path=input_topol_path,
            output_dir=temp_output_dir,
            output_name="step",
            verbose=verbose,
        )

        print("!!!!!!!!!!!1")
        print(mdp_template_path)
        # Run MDrun
        mdrun_outputs = self.mdrun.run(
            input_tpr_path=grompp_output,
            output_name=output_prefix,
            verbose=verbose,
        )

        # Move selected files to main output directory
        files_to_keep = [
            mdrun_outputs[ext] for ext in keep_files if ext in mdrun_outputs
        ]
        output_files = batch_copy_file(
            files_to_keep, main_output_dir, delete_original=True
        )

        logger.info(f"Workflow step completed. Output files saved to {main_output_dir}")
        return output_files
