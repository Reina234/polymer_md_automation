from typing import Dict, List
from A_modules.atomistic.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.shared.utils.file_utils import (
    directory_exists_check_wrapper,
    copy_file,
    rename_file,
)


class FullEquilibrationWorkflow:
    def __init__(self, mdp_cache: MDPCache):
        self.mdp_cache = mdp_cache
        self.steps = []

    def add_step(
        self,
        step_name: str,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
    ):
        self.steps.append((step_name, workflow_step, template_path, base_params))

    @directory_exists_check_wrapper(dir_arg_index=3)
    @directory_exists_check_wrapper(dir_arg_index=4)
    @directory_exists_check_wrapper(dir_arg_index=5)
    def run(
        self,
        input_gro_path: str,
        input_topol_path: str,
        temp_output_dir: str,
        main_output_dir: str,
        log_dir: str,
        final_gro_name: str,
        varying_params_list: List[Dict[str, str]],
        save_intermediate_edr: bool = False,
        save_intermediate_gro: bool = False,
        save_intermediate_log: bool = False,
        verbose: bool = False,
    ):
        current_gro_path = input_gro_path  # Initialize with the initial input GRO file

        for varying_params in varying_params_list:
            for step_name, step, template_path, base_params in self.steps:
                # Merge base and varying parameters
                params = {**base_params, **varying_params}

                print("!!!!!!!!!!!!!!!!!!")
                print(current_gro_path)
                # Run the step
                current_gro_path = step.run(
                    step_name=step_name,
                    mdp_template_path=template_path,
                    input_gro_path=current_gro_path,
                    input_topol_path=input_topol_path,
                    temp_output_dir=temp_output_dir,
                    log_dir=log_dir,
                    varying_params=params,
                    mdp_cache=self.mdp_cache,
                    save_intermediate_edr=save_intermediate_edr,
                    save_intermediate_gro=save_intermediate_gro,
                    save_intermediate_log=save_intermediate_log,
                    verbose=verbose,
                )

                if not current_gro_path:
                    raise RuntimeError(
                        f"Step '{step_name}' did not generate a .gro file required for subsequent steps."
                    )

        final_gro_path = copy_file(
            current_gro_path, main_output_dir, delete_original=False
        )
        final_gro_path = rename_file(final_gro_path, final_gro_name)

        return final_gro_path
