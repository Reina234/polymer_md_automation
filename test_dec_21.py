from A_modules.shared.packmol.solvent_box import PackmolSolventBox
from data_models.solvent import Solvent
import os
from A_modules.shared.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)

gro_file = "TEST/POLY_GMX.gro"

converter = EditconfGROtoPDBConverter()
box_size_nm = [10, 10, 10]
# may need to fix file extensions on converters
output_pdb = converter.run(
    gro_file, "ZZZZZ", "test_output.pdb", box_size_nm=box_size_nm
)


# Path to a valid solvent PDB file
output_dir = "test_output"
output_file = "test_solvent_box.pdb"
# Example box size in nanometers
solvent = Solvent("Hexane", 186.18, 660, output_pdb, "TMZK")

# Create an instance of PackmolSolventBox
packmol_operation = PackmolSolventBox()

# Run the Packmol operation
result = packmol_operation.run(
    output_pdb,
    output_file,
    output_dir="ZZZZZ",
    solvent=solvent,
    box_size_nm=box_size_nm,
)

# Validate the result
assert os.path.exists(result), f"Output file not found: {result}"


print(result)
from A_modules.shared.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)

converter = EditconfPDBtoGROConverter()
output_gro = converter.run(result, "ZZZZZ", "test_output.gro")


from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.gromacs.equilibriation.workflow_step.npt_workflow_step import (
    NptWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun

# Initialize MDP cache
mdp_cache = MDPCache(cache_dir="cache")

# Set up workflow
workflow = FullEquilibrationWorkflow(mdp_cache)

workflow.add_step(
    workflow_step=NptWorkflowStep(Grompp(), MDrun()),
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/em.mdp",
    base_params=None,
)

# List of varying parameters
# varying_params_list = [{"temp": str(temp), "pressure_tau": "2.0"} for temp in [300]]

# Run workflow
workflow.run(
    input_gro_path=output_gro,
    input_topol_path="/temp/TMZK.acpype/TMZK_GMX.top",
    temp_output_dir="temp_outputs",
    main_output_dir="final_outputs",
    keep_files=["gro", "log"],
    varying_params_list=None,
    verbose=True,
)


# NEED TO ADD BOX SIZE INTO EDITCONF SOMEHOW!
