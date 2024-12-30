from typing import Dict, List
from modules.atomistic.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from modules.shared.utils.file_utils import (
    directory_exists_check_wrapper,
    copy_file,
    rename_file,
)
import os
from modules.atomistic.utils.mdp_utils import generate_dynamic_filename


class FullEquilibrationWorkflow:
    def __init__(self, mdp_cache: MDPCache):
        self.mdp_cache = mdp_cache
        self.em_steps = []  # Store EM steps separately
        self.thermal_steps = []  # Store temperature-dependent steps

    def add_em_step(
        self,
        step_name: str,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
    ):
        """Add energy minimization step."""
        self.em_steps.append((step_name, workflow_step, template_path, base_params))

    def add_thermal_step(
        self,
        step_name: str,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
    ):
        """Add thermal steps (e.g., NVT, NPT)."""
        self.thermal_steps.append(
            (step_name, workflow_step, template_path, base_params)
        )

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
        varying_params_list: List[Dict[str, str]],
        override_safeguard_off: bool = False,
        save_intermediate_edr: bool = False,
        save_intermediate_gro: bool = False,
        save_intermediate_log: bool = False,
        verbose: bool = False,
    ):
        subdir = "equilibriated_gros"
        main_output_dir = os.path.join(main_output_dir, subdir)
        current_gro_path = input_gro_path

        # Run all EM steps first
        for step_name, step, template_path, base_params in self.em_steps:
            current_gro_path = step.run(
                step_name=step_name,
                mdp_template_path=template_path,
                input_gro_path=current_gro_path,
                input_topol_path=input_topol_path,
                temp_output_dir=temp_output_dir,
                log_dir=log_dir,
                varying_params=base_params,  # No variation for EM steps
                mdp_cache=self.mdp_cache,
                save_intermediate_edr=save_intermediate_edr,
                save_intermediate_gro=save_intermediate_gro,
                save_intermediate_log=save_intermediate_log,
                verbose=verbose,
            )

            if not current_gro_path:
                raise RuntimeError(
                    f"EM Step '{step_name}' did not generate a .gro file required for subsequent steps."
                )

        # Run thermal steps with varying parameters
        for varying_params in varying_params_list:
            for step_name, step, template_path, base_params in self.thermal_steps:
                # Merge base and varying parameters
                params = {**base_params, **varying_params}

                # Run the step
                current_gro_path = step.run(
                    step_name=f"{step_name}",
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

            final_gro_name = generate_dynamic_filename(varying_params)
            final_gro_path = copy_file(
                current_gro_path,
                main_output_dir,
                delete_original=override_safeguard_off,
            )
            rename_file(final_gro_path, final_gro_name)

        return main_output_dir, final_gro_path
