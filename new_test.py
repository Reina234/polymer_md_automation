from modules.atomistic.workflows.parameterised_polymer_generator.full_workflow import (
    run_polymerisation_workflow,
)

from modules.atomistic.workflows.equilibriated_solvent_box_generator.full_workflow import (
    run_solvent_workflow,
)
from modules.atomistic.workflows.solvated_neutralised_polymer.full_workflow import (
    run_polymer_solvation_workflow,
)
from data_models.solvent import Solvent

polymer_paths = run_polymerisation_workflow("C=Cc1ccccc1", 3, "1_22_test")
solvent = Solvent("hexane", 86.17, 661, "input_data/solvent_pdbs/hexane.pdb", 1.24e-4)

outputs = run_solvent_workflow(
    solvent=solvent,
    box_size_nm=[5, 5, 5],
    temperatures=[298],
    output_dir="1_22_test",
)

run_polymer_solvation_workflow(
    polymer_paths,
    outputs.gro,
    outputs.itp,
    "1_22_test",
    "box",
    298,
    solvent.compressibility,
)
