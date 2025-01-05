from modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)

final_dir, outputs = workflow.run(
    "temp/polymer_in_solvent.gro",
    "temp/topol.top",
    "temp",
    "test",
    "logs",
    varying_params_list=[{"temp": str(298), "compressibility": 1.24e-4}],
    files_to_keep=["edr", "trr", "gro", "xtc"],
    subdir="gromacs_output",
    verbose=False,
)
# keep topol file
