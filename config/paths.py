import os
from enum import Enum

BASE_OUTPUT_DIR = "output"
ACPYPE_PARAMETERIZER_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "acpype_output")
ACPYPE_BASE_NAME = "POLY"
GROMACS_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "gromacs_output")


TEMPORARY_OUTPUT_DIR = "temp"

ACPYPE_PARAMETERIZER_OUTPUT_SUBDIR = "acpype_output"
GROMACS_OUTPUT_SUBDIR = "gromacs"
TOPOL_NAME = "topol.top"


MDP_DIRS = {
    "nvt": "gromacs/mdp_scripts/nvt",
    "npt": "gromacs/mdp_scripts/npt",
}

MDP_NAMING_SCHEME = "{run_type}_{temp}K.mdp"
MDP_TEMPLATE_PATHS = {
    "nvt": "input/templates/gromacs_mdp/nvt.mdp",
    "npt": "input/templates/gromacs_mdp/npt.mdp",
}


class TemplatedMdps(Enum):
    NVT = "nvt"
    NPT = "npt"
