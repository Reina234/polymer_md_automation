from A_modules.atomistic.gromacs.gromacs_config import MDP_NAMING_SCHEME, TemplatedMdps
from A_modules.shared.utils.utils import (
    directory_exists_check_wrapper,
    generate_file_from_template,
)
from typing import Optional
import os
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@directory_exists_check_wrapper(dir_arg_index=1)
def get_or_create_mdp_file(
    mdp_template: TemplatedMdps,
    output_dir: str,
    simulation_temp_k: float,
    mdp_output_naming_scheme: str = MDP_NAMING_SCHEME,
):
    template_path = mdp_template.template_path
    output_filename = mdp_template.generate_mdp_name(
        temp=simulation_temp_k, mdp_output_naming_scheme=mdp_output_naming_scheme
    )
    output_path = os.path.join(output_dir, output_filename)
    generate_file_from_template(
        template_path,
        output_path,
        replacements={"temp": str(simulation_temp_k)},
        overwrite=False,
    )
    return output_path
