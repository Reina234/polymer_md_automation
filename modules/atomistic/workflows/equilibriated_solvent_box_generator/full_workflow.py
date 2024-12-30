from data_models.solvent import Solvent
from modules.atomistic.utils.file_utils import create_solvent_box_gro
from modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from modules.shared.file_conversion.converter_factory import ConverterFactory
from modules.atomistic.workflows.equilibriated_solvent_box_generator.file_preparation_utils import (
    process_solvent_files,
)
from modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from config.paths import EQUILIBRIATED_SOLVENT_BOX_DIR
from typing import List, Optional
from data_models.output_types import GromacsPaths
from config.paths import TEMP_DIR, LOG_DIR
from modules.shared.utils.file_utils import add_identifier_name, delete_directory
from modules.atomistic.utils.mdp_utils import format_temperatures

import os


def run_solvent_workflow(
    solvent: Solvent,
    box_size_nm: List[float],
    workflow: FullEquilibrationWorkflow,
    temperatures: List[float],
    output_dir: str = EQUILIBRIATED_SOLVENT_BOX_DIR,
    identifier: Optional[str] = None,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
    override_safeguard_off=False,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    cleanup: bool = True,
    confirm_log_deletion: bool = True,
) -> str:
    subdir = add_identifier_name(solvent.name, identifier=identifier, suffix=None)
    output_dir = os.path.join(output_dir, subdir)
    parameterizer = parameterizer(acpype_molecule_name=solvent.pdb_molecule_name)
    itp_name = "solvent"
    mol2_converter = ConverterFactory().get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(solvent.pdb_path, temp_dir, verbose=verbose)
    file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
    parametized_files = parameterizer.run(
        mol2_file, temp_dir, file_config, verbose=verbose
    )
    solvent_box_gro = create_solvent_box_gro(
        parametized_files.gro_path, temp_dir, box_size_nm=box_size_nm, solvent=solvent
    )
    reformatted_files = process_solvent_files(
        parametized_files.itp_path,
        solvent_box_gro,
        parametized_files.top_path,
        new_residue_name=solvent.pdb_molecule_name,
        output_itp_dir=output_dir,
        output_gro_dir=temp_dir,
        output_topol_dir=temp_dir,
        output_itp_name=itp_name,
    )

    temperatures = format_temperatures(temperatures)

    gro_name = add_identifier_name(
        solvent.name,
        default_identifier="equilibriated",
        identifier=identifier,
        suffix=None,
    )
    gro_output_dir, _ = workflow.run(
        reformatted_files.gro_path,
        reformatted_files.top_path,
        temp_dir,
        output_dir,
        log_dir,
        override_safeguard_off=override_safeguard_off,
        varying_params_list=temperatures,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        verbose=verbose,
    )

    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=False)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)

    return gro_output_dir
