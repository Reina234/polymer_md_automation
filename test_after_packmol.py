from A_modules.shared.packmol.solvent_box import PackmolSolventBox
from data_models.solvent import Solvent
import os
from A_modules.shared.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)


from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.gromacs.equilibriation.workflow_step.npt_workflow_step import (
    NptWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun

output_gro = "ZZZZZ/test2_solvent_box.gro"
# Initialize MDP cache
mdp_cache = MDPCache(cache_dir="cache")

# Set up workflow
workflow = FullEquilibrationWorkflow(mdp_cache)

workflow.add_step(
    workflow_step=NptWorkflowStep(Grompp(), MDrun()),
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/npt.mdp",
    base_params={"pressure": "1.0", "compressibility": "4.5e-5"},
)

# List of varying parameters
varying_params_list = [{"temp": str(temp), "pressure_tau": "2.0"} for temp in [300]]

# Run workflow
workflow.run(
    input_gro_path=output_gro,
    input_topol_path="temp/TMZK.acpype/TMZK_GMX.top",
    temp_output_dir="temp_outputs",
    main_output_dir="final_outputs",
    keep_files=["gro", "log"],
    varying_params_list=varying_params_list,
    verbose=True,
)


# NEED TO ADD BOX SIZE INTO EDITCONF SOMEHOW!
