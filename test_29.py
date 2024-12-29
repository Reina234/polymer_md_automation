from solute_preparation_utils import prepare_solute_files

solvent_itp_file = "TEST_RUN_28/hexane_run_1/solvent.itp"
solute_itp_file = "temp_29_test/POLY_GMX.itp"
top_file = "temp_29_test/POLY_GMX.top"
solvent_box_gro = "temp_29_test/solute_box.gro"
result = prepare_solute_files(
    solute_itp_file=solute_itp_file,
    solvent_itp_file=solvent_itp_file,
    solvent_box_gro_file=solvent_box_gro,
    input_top_file=top_file,
    output_dir="temp_29_test/prepared_files",
)

gro = result.gro_path
itp = result.itp_path
top = result.top_path
from A_modules.atomistic.gromacs.workflows.equilibriated_solvent_box.workflow_config import (
    workflow,
)
