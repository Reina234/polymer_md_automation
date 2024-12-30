from A_modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)

workflow.run(
    "temp_29_test/prepared_files/solute_box.gro",
    "temp_29_test/prepared_files/topol.top",
    "temp_29_test/temp",
    "temp_29_test/equilibriation_run",
    "logs_29",
    [{"temp": "298"}],
)
