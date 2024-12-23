from typing import Dict, List
from A_modules.atomistic.gromacs.equilibriation.workflow_step.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache


class FullEquilibrationWorkflow:
    def __init__(self, mdp_cache: MDPCache):
        self.mdp_cache = mdp_cache
        self.steps = []

    def add_step(
        self,
        workflow_step: BaseWorkflowStep,
        template_path: str,
        base_params: Dict[str, str],
    ):
        self.steps.append((workflow_step, template_path, base_params))

    def run(
        self,
        input_gro_path: str,
        input_topol_path: str,
        temp_output_dir: str,
        main_output_dir: str,
        keep_files: List[str],
        varying_params_list: List[Dict[str, str]],
        verbose: bool = False,
    ):
        for varying_params in varying_params_list:
            current_gro_path = (
                input_gro_path  # Initialize with the initial input GRO file
            )

            for step, template_path, base_params in self.steps:
                # Merge base and varying parameters
                params = {**base_params, **varying_params}

                # Run the step
                output_files = step.run(
                    mdp_template_path=template_path,
                    input_gro_path=current_gro_path,
                    input_topol_path=input_topol_path,
                    temp_output_dir=temp_output_dir,
                    main_output_dir=main_output_dir,
                    keep_files=keep_files,
                    varying_params=params,
                    mdp_cache=self.mdp_cache,
                    verbose=verbose,
                )

                # Update the current GRO file path for the next step
                current_gro_path = next(
                    (file for file in output_files if file.endswith(".gro")), None
                )
                print("!!!!!!!!!!!")
                print(f"Current GRO path: {current_gro_path}")
                if not current_gro_path:
                    raise RuntimeError(
                        f"Step did not generate a .gro file required for subsequent steps. Check the step configuration."
                    )
