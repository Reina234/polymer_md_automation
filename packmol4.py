from A_modules.shared.packmol.solvent_box import PackmolSolventBox
from data_models.solvent import Solvent
import os


def test_packmol_run():
    # Define test parameters
    solvent_file = "input/solvents/pdb/hexane.pdb"  # Path to a valid solvent PDB file
    output_dir = "test_output"
    output_file = "test_solvent_box.pdb"
    box_size_nm = [10, 10, 10]  # Example box size in nanometers
    solvent = Solvent("Hexane", 186.18, 660, solvent_file, "TMZK")

    # Create an instance of PackmolSolventBox
    packmol_operation = PackmolSolventBox()

    # Run the Packmol operation
    result = packmol_operation.run(
        solvent_file,
        output_file,
        output_dir="ZZZ",
        solvent=solvent,
        box_size_nm=box_size_nm,
    )

    # Validate the result
    assert os.path.exists(result), f"Output file not found: {result}"
    return result


result = test_packmol_run()
print(result)
from A_modules.shared.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)

converter = EditconfPDBtoGROConverter()
output_gro = converter.run(result, "ZZZ", "test_output.gro")


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
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/npt.mdp",
    base_params={"pressure": "1.0", "compressibility": "4.5e-5"},
)

# List of varying parameters
varying_params_list = [{"temp": str(temp), "pressure_tau": "2.0"} for temp in [300]]

# Run workflow
workflow.run(
    input_gro_path=output_gro,
    input_topol_path="TRIAL/acpype_output/POLY_GMX.top",
    temp_output_dir="temp_outputs",
    main_output_dir="final_outputs",
    keep_files=["gro", "log"],
    varying_params_list=varying_params_list,
    verbose=True,
)
