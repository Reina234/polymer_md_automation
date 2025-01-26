from modules.atomistic.workflows.parameterised_polymer_generator.full_workflow import (
    run_polymerisation_workflow,
)

from modules.atomistic.workflows.equilibriated_solvent_box_generator.full_workflow import (
    run_solvent_workflow,
)
from modules.atomistic.workflows.solvated_neutralised_polymer.workflow_no_mdrun import (
    run_polymer_solvation_workflow,
)
from data_models.solvent import Solvent
from data_models.output_types import GromacsPaths

polymer_paths = GromacsPaths(
    itp_path="test_3.itp",
    gro_path="manual_gro.gro",
    top_path="topol_manual.top",
)
solvent = Solvent("hexane", 86.17, 661, "input_data/solvent_pdbs/hexane.pdb", 1.24e-4)


run_polymer_solvation_workflow(
    polymer_paths,
    "1_22_test/hexane/equilibriated_outputs/temp_298_compressibility_0.000124.gro",
    "1_22_test/hexane/solvent.itp",
    "1_22_test",
    "box",
    298,
    solvent.compressibility,
)
