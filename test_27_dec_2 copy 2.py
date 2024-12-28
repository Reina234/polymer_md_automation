from data_models.output_types import GromacsPaths

from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.gromacs.equilibriation.base_workflow_step import (
    BaseWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun

output_files = GromacsPaths(
    "Z2/POLY_GMX.itp",
    "TEST_RUN_27_2/hexane_solvent_box.gro",
    "TEST_RUN_27/POLY_GMX.top",
    None,
)


mdp_cache = MDPCache(cache_dir="cache")
workflow_step = BaseWorkflowStep(Grompp(), MDrun())
workflow = FullEquilibrationWorkflow(mdp_cache)

# Add steps
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

varying_params_list = [{"temp": str(temp)} for temp in [300]]

# Run workflow
workflow.run(
    output_files.gro_path,
    output_files.top_path,
    "temp",
    "output_27_dec",
    "log_dir",
    varying_params_list=varying_params_list,  # Example with no varying parameters
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    verbose=True,
)
