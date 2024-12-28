from A_config.paths import CACHE_DIR
from A_modules.atomistic.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache

mdp_cache = MDPCache(cache_dir=CACHE_DIR)
workflow_step = BaseWorkflowStep(Grompp(), MDrun())
workflow = FullEquilibrationWorkflow(mdp_cache)


workflow.add_step(
    step_name="minim",
    workflow_step=workflow_step,
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/em.mdp",
    base_params={"nsteps": "50000"},
)
workflow.add_step(
    step_name="npt",
    workflow_step=workflow_step,
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/npt.mdp",
    base_params={"pressure": "1.0", "compressibility": "4.5e-5", "pressure_tau": "2.0"},
)
