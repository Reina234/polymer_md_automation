from modules.atomistic.workflows.parameterised_polymer_generator.full_workflow import (
    run_polymerisation_workflow,
)
from modules.atomistic.workflows.solvated_polymer_generator.full_workflow import (
    run_polymer_solvation_workflow,
)
from modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)
from data_models.output_types import GromacsPaths


# mol2_file = run_polymerisation_workflow("C=Cc1ccccc1", 5, "rdkit_test")
polymerised_files = GromacsPaths(
    itp_path="rdkit_test/POLY_GMX.itp",
    gro_path="rdkit_test/POLY_GMX.gro",
    top_path="rdkit_test/POLY_GMX.top",
)

solvent_itp_file = "TEST_RUN_28/hexane_run_1/solvent.itp"
solvent_box_gro = "TEST_RUN_28/hexane_run_1/equilibriated_gros/temp_298.gro"
run_polymer_solvation_workflow(
    polymerised_files,
    subdir="styrene",
    solvent_equilibriated_gro=solvent_box_gro,
    solvent_itp=solvent_itp_file,
    output_dir="test_29_full",
    workflow=workflow,
    temperature=298,
)

############################################################################################
# need to pass in compressibility somehow, set a default, put into Solvent class
# MAYBE, use the file name itself to get the compressibility value and the temp value? !!!
############################################################################################
