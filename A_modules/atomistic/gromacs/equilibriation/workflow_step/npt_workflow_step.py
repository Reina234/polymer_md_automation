from A_modules.atomistic.gromacs.equilibriation.workflow_step.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.shared.utils.file_utils import check_directory_exists, batch_copy_file
from typing import Dict, List, Optional
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
import os


from A_modules.atomistic.gromacs.equilibriation.workflow_step.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from typing import Dict, Optional


class NptWorkflowStep(BaseWorkflowStep):
    def __init__(
        self,
        grompp: Grompp,
        mdrun: MDrun,
        mdp_cache: Optional[MDPCache] = None,
        default_params: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize the NPT Workflow Step.

        :param grompp: Grompp command instance.
        :param mdrun: MDrun command instance.
        :param mdp_cache: MDP cache for managing and reusing MDP files.
        :param default_params: Default parameters specific to NPT runs.
        """
        super().__init__(
            grompp,
            mdrun,
            default_params
            or {
                "pressure": "1.0",
                "compressibility": "4.5e-5",
            },
        )

    def merge_params(self, varying_params: Dict[str, str]) -> Dict[str, str]:
        """
        Merge varying parameters with default NPT parameters.
        Priority is given to varying parameters.

        :param varying_params: Parameters to override the defaults.
        :return: Merged parameters.
        """
        return super().merge_params(varying_params)
