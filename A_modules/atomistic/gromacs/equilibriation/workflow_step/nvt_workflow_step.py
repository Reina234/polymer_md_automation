from A_modules.atomistic.gromacs.equilibriation.workflow_step.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.shared.utils.file_utils import check_directory_exists, batch_copy_file
from typing import Dict, List
import os


class NvtWorkflowStep(BaseWorkflowStep):
    def run(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        temp_output_dir: str,
        main_output_dir: str,
        keep_files: List[str],
        verbose: bool = False,
    ) -> Dict[str, str]:
        # Ensure directories exist
        temp_output_dir = check_directory_exists(temp_output_dir)
        main_output_dir = check_directory_exists(main_output_dir)

        # Run GROMPP
        tpr_file = self.grompp.run(
            mdp_file_path=mdp_file_path,
            input_gro_path=input_gro_path,
            input_topol_path=input_topol_path,
            output_tpr_dir=temp_output_dir,
            output_tpr_name="nvt",
            verbose=verbose,
        )

        # Run MDrun
        mdrun_outputs = self.mdrun.run(
            input_tpr_path=tpr_file,
            output_name=os.path.join(temp_output_dir, "nvt"),
            verbose=verbose,
        )

        # Move kept files
        kept_files = {
            ext: mdrun_outputs[ext]
            for ext in keep_files
            if ext in mdrun_outputs and mdrun_outputs[ext]
        }
        batch_copy_file(
            files_to_move=list(kept_files.values()),
            dest_dir=main_output_dir,
            delete_original=True,
        )

        return {
            ext: os.path.join(main_output_dir, os.path.basename(path))
            for ext, path in kept_files.items()
        }
