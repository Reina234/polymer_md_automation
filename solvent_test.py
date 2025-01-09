from modules.atomistic.workflows.equilibriated_solvent_box_generator.full_workflow import (
    run_solvent_workflow,
)
from modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)
from data_models.solvent import Solvent

solvent = Solvent("hexane", 86.17, 661, "input_data/solvent_pdbs/hexane.pdb", 1.24e-4)

run_solvent_workflow(
    solvent=solvent,
    box_size_nm=[5, 5, 5],
    temperatures=[298, 318],
    output_dir="1_5_test",
)
