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
workflow2 = FullEquilibrationWorkflow(mdp_cache)


workflow2.add_em_step(
    step_name="minim_1",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/em_2.mdp",
    base_params={"nsteps": "50000", "emtol": "1000"},
)
workflow2.add_em_step(
    step_name="minim_2",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/em_2.mdp",
    base_params={"nsteps": "50000", "emtol": "100"},
)

workflow2.add_thermal_step(
    step_name="nvt_short",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/nvt.mdp",
    base_params={
        "nsteps": "100000",
    },
)

workflow2.add_thermal_step(
    step_name="npt_short_1",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/npt_short.mdp",
    base_params={
        "nsteps": "50000",
        "dt": "0.001",
    },
)
workflow2.add_thermal_step(
    step_name="npt_short_2",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/npt_short.mdp",
    base_params={
        "nsteps": "30000",
        "dt": "0.002",
    },
)
workflow2.add_thermal_step(
    step_name="npt_2",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/npt_2.mdp",
    base_params={
        "nsteps": "10000",
        "dt": "0.001",
    },
)

workflow2.add_thermal_step(
    step_name="production",
    workflow_step=workflow_step,
    template_path="modules/atomistic/gromacs/equilibriation/templates/ion_mdps/prod.mdp",
    base_params={
        "nsteps": "10000",
        "dt": "0.002",
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
