import os
from enum import Enum

BASE_OUTPUT_DIR = "output"
ACPYPE_PARAMETERIZER_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "acpype_output")
ACPYPE_BASE_NAME = "POLY"
GROMACS_OUTPUT_DIR = os.path.join(BASE_OUTPUT_DIR, "gromacs_output")


class FilesToExtract(Enum):
    GRO = "gro"
    TOP = "top"
    ITP = "itp"


TEMPORARY_OUTPUT_DIR = "temp"

ACPYPE_PARAMETERIZER_OUTPUT_SUBDIR = "acpype_output"
GROMACS_OUTPUT_SUBDIR = "gromacs"
TOPOL_NAME = "topol.top"
EQUILIBRIUM_SUBDIR = "equilibrium"
# saves to output/runname/gromacs/equilibrium


SOLVENT_ITP_DIR = "preprocessed_output/solvents/solvent_itp"
SOLVENT_MOL2_DIR = "preprocessed_output/solvents/solvent_mol2"
SOLVENT_JSON_PATH = "preprocessed_output/solvents/solvent_atomtypes.json"

MDP_DIRS = {
    "nvt": "preprocessed_output/gromacs/nvt",
    "npt": "preprocessed_output/gromacs/npt",
}

MDP_NAMING_SCHEME = "{run_type}_{temp}K.mdp"
MDP_TEMPLATE_PATHS = {
    "nvt": "input/templates/gromacs_mdp/nvt.mdp",
    "npt": "input/templates/gromacs_mdp/npt.mdp",
}

# NOTE!! THIS IS DEPENDENT ON THE MDP NAMING SCHEME AND MDP DIRS!!########
MDP_FULL_PATHS = {
    "nvt": "preprocessed_output/gromacs/nvt/nvt_{temp}K.mdp",
    "npt": "preprocessed_output/gromacs/npt/npt_{temp}K.mdp",
    "minim": "preprocessed_output/gromacs/minim.mdp",
}
###########################################################################


class TemplatedMdps(Enum):
    NVT = "nvt"
    NPT = "npt"
    MINIM = "minim"


ACPYPE_SOLVENT_OUTPUT_SUBDIR = "acpype_solvent_output"
