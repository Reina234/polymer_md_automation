from A_modules.atomistic.workflows.parameterised_polymer_generator.full_workflow import (
    run_polymerisation_workflow,
)
from A_modules.atomistic.workflows.solvated_polymer_generator.full_workflow import (
    run_polymer_solvation_workflow,
)
from A_modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)
from data_models.output_types import GromacsPaths


# mol2_file = run_polymerisation_workflow("C=Cc1ccccc1", 5, "rdkit_test")
polymerised_files = GromacsPaths(
    itp_path="rdkit_test/POLY_GMX.itp",
    gro_path="trial.gro",
    top_path="rdkit_test/POLY_GMX.top",
)

solvent_itp_file = "TEST_RUN_28/hexane_run_1/solvent.itp"
solvent_box_gro = "temp_29_test/solute_box.gro"
run_polymer_solvation_workflow(
    polymerised_files,
    solvent_equilibriated_gro=solvent_box_gro,
    solvent_itp=solvent_itp_file,
    output_dir="test_29_full",
    workflow=workflow,
    temperature=298,
)
