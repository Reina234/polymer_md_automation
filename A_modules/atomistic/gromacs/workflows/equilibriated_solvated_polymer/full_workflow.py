from typing import List, Optional
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from A_config.paths import TEMP_DIR, LOG_DIR
from A_modules.atomistic.gromacs.commands.insert_molecules import InsertMolecules
from A_modules.atomistic.gromacs.workflows.equilibriated_solvated_polymer.file_preparation_utils import (
    prepare_solute_files,
)
from data_models.output_types import GromacsPaths


# NOTE: will need helper function to extract temperatures, solvent gro, solvent itp, all from folder
def run_polymer_solvation_workflow(
    parameterised_polymer: GromacsPaths,
    solvent_equilibriated_gro: str,
    solvent_itp: str,
    output_dir: str,
    workflow: FullEquilibrationWorkflow,
    temperature: float,
    polymer_name: str = "POLY",
    num_polymers: int = 1,
    identifier: Optional[str] = None,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
    override_safeguard_off=False,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
):

    polymer_in_solvent = InsertMolecules().run(
        solvent_equilibriated_gro,
        parameterised_polymer.gro_path,
        num_polymers,
        temp_dir,
        "polymer_in_solvent",
    )

    prepared_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=polymer_in_solvent,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    final_gro = workflow.run(
        prepared_files.gro_path,
        prepared_files.top_path,
        temp_output_dir=temp_dir,
        main_output_dir=output_dir,
        log_dir=log_dir,
        varying_params_list=[{"temp": temp} for temp in temperature],
        override_safeguard_off=override_safeguard_off,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        verbose=verbose,
    )

    return final_gro
