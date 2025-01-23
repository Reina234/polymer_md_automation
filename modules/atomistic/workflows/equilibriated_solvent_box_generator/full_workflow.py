from data_models.solvent import Solvent
from config.acpype_config import AcpypeOutputConfig
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
from modules.atomistic.workflows.equilibriated_solvent_box_generator.workflow_config import (
    workflow,
)
from modules.atomistic.utils.file_utils import (
    get_gro_handler,
    get_residue_number,
    rename_residue_name_from_handler,
    create_solvent_box_gro,
    rename_data_column_content,
    export_gro_handler,
    create_includes_section,
    delete_all_include_sections,
)
from data_models.output_types import GromacsPaths
from config.paths import EQUILIBRIATED_SOLVENT_BOX_DIR
from typing import List, Optional
from data_models.output_types import GromacsOutputs
from config.paths import TEMP_DIR, LOG_DIR
from modules.shared.utils.file_utils import add_identifier_name, delete_directory
from modules.atomistic.utils.mdp_utils import format_temperatures
from modules.atomistic.workflows.equilibriated_solvent_box_generator.full_workflow import (
    workflow,
)


import os


def run_solvent_workflow(
    solvent: Solvent,
    box_size_nm: List[float],
    temperatures: List[float],
    output_dir: str = EQUILIBRIATED_SOLVENT_BOX_DIR,
    identifier: Optional[str] = None,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    cleanup: bool = True,
    confirm_log_deletion: bool = True,
) -> GromacsOutputs:

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

    varying_params_list = format_temperatures(
        temperatures=temperatures, compressibility=solvent.compressibility
    )

    gro_name = add_identifier_name(
        solvent.name,
        default_identifier="equilibriated",
        identifier=identifier,
        suffix=None,
    )
    gro_output_dir, output_paths = workflow.run(
        reformatted_files.gro_path,
        reformatted_files.top_path,
        temp_dir,
        output_dir,
        log_dir,
        files_to_keep=["gro", "tpr"],
        varying_params_list=varying_params_list,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        verbose=verbose,
    )

    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=False)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)

    output_paths.itp = reformatted_files.itp_path
    return output_paths
