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

# NOTE!! THIS IS DEPENDENT ON THE MDP NAMING SCHEME AND MDP DIRS!!########
MDP_FULL_PATHS = {
    "nvt": "gromacs/mdp_scripts/nvt/nvt_{temp}K.mdp",
    "npt": "gromacs/mdp_scripts/npt/npt_{temp}K.mdp",
    "minim": "gromacs/mdp_scripts/minim.mdp",
}
###########################################################################


class TemplatedMdps(Enum):
    NVT = "nvt"
    NPT = "npt"
    MINIM = "minim"


ACPYPE_SOLVENT_OUTPUT_SUBDIR = "acpype_solvent_output"
