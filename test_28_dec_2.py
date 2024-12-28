from A_modules.atomistic.gromacs.workflows.equilibriated_solvent_box.workflow_config import (
    workflow,
)
from A_modules.atomistic.gromacs.workflows.equilibriated_solvent_box.full_workflow import (
    run_solvent_workflow,
)
from data_models.solvent import Solvent

solvent = Solvent("Hexane", 86, 660, "input/solvents/pdb/hexane.pdb", "TMZK")
box_size_nm = [5, 5, 5]

result = run_solvent_workflow(
    solvent=solvent,
    box_size_nm=box_size_nm,
    output_dir="TEST_RUN_28",
    workflow=workflow,
    temperatures=[298, 318],
    identifier="run_1",
    verbose=False,
)
