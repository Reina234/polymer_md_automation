from config.paths import MDP_DIRS, MDP_NAMING_SCHEME, MDP_TEMPLATE_PATHS, TemplatedMdps
from typing import Optional
import os
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def create_mdps(
    mdp_type: TemplatedMdps,
    simulation_temp_k: float,
    template_path: Optional[str] = None,
    output_dir: Optional[str] = None,
):
    if not template_path:
        template_path = MDP_TEMPLATE_PATHS[mdp_type.value]

    if not output_dir:
        output_dir = MDP_DIRS[mdp_type.value]

    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template path {template_path} does not exist")

    with open(template_path, "r") as file:
        content = file.read()

    content = content.replace("{temp}", str(simulation_temp_k))
    output_filename = MDP_NAMING_SCHEME.format(
        run_type=mdp_type.value, temp=simulation_temp_k
    )
    output_path = os.path.join(output_dir, output_filename)
    os.makedirs(output_dir, exist_ok=True)

    with open(output_path, "w") as file:
        file.write(content)

    logger.info(f"Created MDP file at {output_path}")
    return output_path
