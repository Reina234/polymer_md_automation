from test_28_dec import workflow
from solvent_workflow import solvent_workflow
from data_models.solvent import Solvent

solvent = Solvent("Hexane", 86, 660, "input/solvents/pdb/hexane.pdb", "TMZK")
box_size_nm = [5, 5, 5]

result = solvent_workflow(
    solvent=solvent,
    box_size_nm=box_size_nm,
    output_dir="TEST_RUN_28",
    workflow=workflow,
    temperatures=[300],
    identifier="run_1",
)
