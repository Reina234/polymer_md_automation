from typing import List, Optional
from modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from config.paths import TEMP_DIR, LOG_DIR
from modules.atomistic.gromacs.commands.insert_molecules import InsertMolecules
from modules.atomistic.workflows.solvated_polymer_generator.file_preparation_utils import (
    prepare_solute_files,
)
from modules.shared.utils.file_utils import delete_directory, copy_file
from data_models.output_types import GromacsPaths
from modules.atomistic.utils.moltemplate_utils import add_polymer_to_solvent
import os


# NOTE: will need naming scheme
# NOTE: will need helper function to extract temperatures, solvent gro, solvent itp, all from folder
def run_polymer_solvation_workflow(
    parameterised_polymer: GromacsPaths,
    solvent_equilibriated_gro: str,
    solvent_itp: str,
    output_dir: str,
    subdir: str,
    workflow: FullEquilibrationWorkflow,
    temperature: float,
    polymer_name: str = "POLY",
    cutoff: float = 0.2,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
    override_safeguard_off=False,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    cleanup: bool = True,
    confirm_log_deletion: bool = True,
) -> str:
    final_output_dir = os.path.join(output_dir, subdir)
    polymer_in_solvent = add_polymer_to_solvent(
        parameterised_polymer.gro_path,
        solvent_equilibriated_gro,
        temp_dir,
        "polymer_in_solvent",
        cutoff=cutoff,
    )

    prepared_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=polymer_in_solvent,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    final_dir, outputs = workflow.run(
        prepared_files.gro_path,
        prepared_files.top_path,
        temp_dir,
        final_output_dir,
        log_dir,
        varying_params_list=[{"temp": str(temperature)}],
        files_to_keep=["edr", "trr", "gro"],
        override_safeguard_off=override_safeguard_off,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        subdir="gromacs_output",
        verbose=verbose,
    )
    # keep topol file
    topol_file = copy_file(
        prepared_files.top_path, final_output_dir, skip_if_exists=True
    )
    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=False)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)
    outputs.top = topol_file
    return outputs


# NOTE: might need more than just gro as final output!
