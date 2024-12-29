from A_modules.atomistic.gromacs.workflows.parameterised_polymer_generator.full_workflow import (
    run_polymerisation_workflow,
)
from A_modules.atomistic.gromacs.workflows.equilibriated_solvated_polymer.full_workflow import (
    run_polymer_solvation_workflow,
)
from A_modules.atomistic.gromacs.workflows.equilibriated_solvent_box.workflow_config import (
    workflow,
)

mol2_file = run_polymerisation_workflow("C=Cc1ccccc1", 5, "rdkit_test")

solvent_itp_file = "TEST_RUN_28/hexane_run_1/solvent.itp"
solvent_box_gro = "temp_29_test/solute_box.gro"
run_polymer_solvation_workflow(
    mol2_file,
    solvent_equilibriated_gro=solvent_box_gro,
    solvent_itp=solvent_itp_file,
    output_dir="test_29_full",
    workflow=workflow,
    temperature=298,
)
