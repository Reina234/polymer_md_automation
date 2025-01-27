from typing import List, Optional
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from config.paths import TEMP_DIR, LOG_DIR
from modules.gromacs.commands.genion import GenIon
from modules.atomistic.workflows.solvated_neutralised_polymer.file_preparation_utils import (
    prepare_solute_files,
)
from modules.utils.shared.file_utils import delete_directory, copy_file
from data_models.output_types import GromacsPaths
from modules.moltemplate.moltemplate_utils import (
    add_polymer_to_solvent,
    add_n_parallel_polymers_to_solvent,
)
import os

from modules.atomistic.workflows.solvated_neutralised_polymer.minim_config import (
    minim_workflow,
)
from modules.atomistic.workflows.solvated_neutralised_polymer.workflow_config import (
    workflow,
)


# NOTE: will need naming scheme
# NOTE: will need helper function to extract temperatures, solvent gro, solvent itp, all from folder
def run_polymer_solvation_workflow(
    parameterised_polymer: GromacsPaths,
    solvent_equilibriated_gro: str,
    solvent_itp: str,
    output_dir: str,
    subdir: str,
    temperature: float,
    compressibility: float,
    polymer_name: str = "POLY",
    cutoff: float = 0.2,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
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

    initial_minim_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=polymer_in_solvent,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    _, outputs = minim_workflow.run(
        initial_minim_files.gro_path,
        initial_minim_files.top_path,
        temp_dir,
        temp_dir,
        log_dir,
        varying_params_list=[None],
        files_to_keep=["gro", "tpr"],
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        subdir="gromacs_output",
        verbose=verbose,
        file_name_override="initial_minim",
    )

    neutralised_gro = GenIon().run(
        outputs.gro, outputs.tpr, initial_minim_files.top_path, temp_dir
    )

    prepared_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=neutralised_gro,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    _, outputs = workflow.run(
        prepared_files.gro_path,
        prepared_files.top_path,
        temp_dir,
        final_output_dir,
        log_dir,
        varying_params_list=[
            {"temp": str(temperature), "compressibility": compressibility}
        ],
        files_to_keep=["edr", "trr", "gro", "xtc", "tpr"],
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        subdir="gromacs_output",
        verbose=verbose,
    )
    # keep topol file
    topol_file = copy_file(
        initial_minim_files.top_path, final_output_dir, skip_if_exists=True
    )
    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=confirm_log_deletion)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)
    outputs.top = topol_file
    return outputs


# NOTE: might need more than just gro as final output!


from typing import List, Optional
from modules.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from config.paths import TEMP_DIR, LOG_DIR
from modules.gromacs.commands.genion import GenIon
from modules.atomistic.workflows.solvated_neutralised_polymer.file_preparation_utils import (
    prepare_solute_files,
)
from modules.utils.shared.file_utils import delete_directory, copy_file
from data_models.output_types import GromacsPaths
from modules.moltemplate.moltemplate_utils import add_polymer_to_solvent
import os

from modules.atomistic.workflows.solvated_neutralised_polymer.minim_config import (
    minim_workflow,
)
from modules.atomistic.workflows.solvated_neutralised_polymer.workflow_config import (
    workflow,
)


# NOTE: will need naming scheme
# NOTE: will need helper function to extract temperatures, solvent gro, solvent itp, all from folder
def run_three_polymer_solvation_workflow(
    parameterised_polymer: GromacsPaths,
    solvent_equilibriated_gro: str,
    solvent_itp: str,
    output_dir: str,
    subdir: str,
    temperature: float,
    compressibility: float,
    polymer_name: str = "POLY",
    cutoff: float = 0.2,
    min_polymer_distance: float = 1.2,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = False,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    cleanup: bool = True,
    confirm_log_deletion: bool = True,
) -> str:
    final_output_dir = os.path.join(output_dir, subdir)
    polymer_in_solvent = add_n_parallel_polymers_to_solvent(
        parameterised_polymer.gro_path,
        solvent_equilibriated_gro,
        2,
        temp_dir,
        "polymer_in_solvent",
        cutoff=cutoff,
        min_distance=min_polymer_distance,
    )

    initial_minim_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=polymer_in_solvent,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    _, outputs = minim_workflow.run(
        initial_minim_files.gro_path,
        initial_minim_files.top_path,
        temp_dir,
        temp_dir,
        log_dir,
        varying_params_list=[None],
        files_to_keep=["gro", "tpr"],
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        subdir="gromacs_output",
        verbose=verbose,
        file_name_override="initial_minim",
    )

    neutralised_gro = GenIon().run(
        outputs.gro, outputs.tpr, initial_minim_files.top_path, temp_dir
    )

    prepared_files = prepare_solute_files(
        solute_itp_file=parameterised_polymer.itp_path,
        solvent_itp_file=solvent_itp,
        solvent_box_gro_file=neutralised_gro,
        input_top_file=parameterised_polymer.top_path,
        output_dir=temp_dir,
        solute_molecule_name=polymer_name,
    )

    _, outputs = workflow.run(
        prepared_files.gro_path,
        prepared_files.top_path,
        temp_dir,
        final_output_dir,
        log_dir,
        varying_params_list=[
            {"temp": str(temperature), "compressibility": compressibility}
        ],
        files_to_keep=["edr", "trr", "gro", "xtc", "tpr"],
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        subdir="gromacs_output",
        verbose=verbose,
    )
    # keep topol file
    topol_file = copy_file(
        initial_minim_files.top_path, final_output_dir, skip_if_exists=True
    )
    if cleanup:
        delete_directory(temp_dir, verbose=verbose, confirm=confirm_log_deletion)
        delete_directory(log_dir, verbose=verbose, confirm=confirm_log_deletion)
    outputs.top = topol_file
    return outputs


# NOTE: might need more than just gro as final output!
