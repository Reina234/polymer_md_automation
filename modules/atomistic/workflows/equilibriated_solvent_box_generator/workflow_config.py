from config.paths import MDP_CACHE_DIR
from modules.atomistic.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from modules.atomistic.gromacs.commands.grompp import Grompp
from modules.atomistic.gromacs.commands.mdrun import MDrun
from modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache

mdp_cache = MDPCache(cache_dir=MDP_CACHE_DIR)
workflow_step = BaseWorkflowStep(Grompp(), MDrun())
workflow = FullEquilibrationWorkflow(mdp_cache)


workflow.add_em_step(
    step_name="soft_em",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/em_2.mdp",
    base_params={"nsteps": "10000", "emtol": "1000"},
)
workflow.add_em_step(
    step_name="minim_1",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/em_2.mdp",
    base_params={"nsteps": "50000", "emtol": "1000"},
)
workflow.add_em_step(
    step_name="minim_2",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/em_2.mdp",
    base_params={"nsteps": "50000", "emtol": "100"},
)

workflow.add_thermal_step(
    step_name="nvt_short",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/nvt.mdp",
    base_params={
        "nsteps": "50000",
    },
)

workflow.add_thermal_step(
    step_name="npt_short",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/npt_short.mdp",
    base_params={
        "nsteps": "20000",
        "dt": "0.002",
    },
)
workflow.add_thermal_step(
    step_name="npt_2",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/npt_2.mdp",
    base_params={
        "nsteps": "15000",
        "dt": "0.001",
    },
)

# workflow.add_step(
#    step_name="npt_full",
#    workflow_step=workflow_step,
#    template_path="A_modules/atomistic/gromacs/equilibriation/templates/npt.mdp",
#    base_params={
#        "nsteps": "5000",
#        "compressibility": "1.4e-4",
#        "dt": "0.001",
#        "pressure": "1.0",
#        "pressure_tau": "2.0",
#    },
# )
