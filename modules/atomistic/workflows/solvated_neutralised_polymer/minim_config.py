from config.paths import MDP_CACHE_DIR
from modules.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from modules.gromacs.commands.grompp import Grompp
from modules.gromacs.commands.mdrun import MDrun
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from modules.gromacs.equilibriation.mdp_cache import MDPCache

mdp_cache = MDPCache(cache_dir=MDP_CACHE_DIR)
workflow_step = BaseWorkflowStep(Grompp(), MDrun())
minim_workflow = FullEquilibrationWorkflow(mdp_cache)


minim_workflow.add_em_step(
    step_name="ion_em",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_em.mdp",
    base_params={"nsteps": "50000", "emtol": "1000"},
)
