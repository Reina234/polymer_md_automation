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
from modules.shared.utils.file_utils import delete_directory
from data_models.output_types import GromacsPaths
from modules.atomistic.utils.moltemplate_utils import add_polymer_to_solvent


# NOTE: will need naming scheme
# NOTE: will need helper function to extract temperatures, solvent gro, solvent itp, all from folder
def run_polymer_solvation_workflow(
    parameterised_polymer: GromacsPaths,
    solvent_equilibriated_gro: str,
    solvent_itp: str,
    output_dir: str,
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

    _, final_gro = workflow.run(
        prepared_files.gro_path,
        prepared_files.top_path,
        temp_dir,
        output_dir,
        log_dir,
        varying_params_list=[{"temp": str(temperature)}],
        override_safeguard_off=override_safeguard_off,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        verbose=verbose,
    )
    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=False)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)

    return final_gro


# NOTE: might need more than just gro as final output!
