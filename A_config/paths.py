import os
from enum import Enum
import dataclasses
from typing import Optional

TEMPORARY_OUTPUT_DIR = "temp"
TOPOL_NAME = "topol.top"


@dataclasses
class AcpypeFiles:
    itp: bool
    gro: bool
    top: bool
    posre: bool

    itp_path: Optional[str]
    gro_path: Optional[str]
    top_path: Optional[str]
    posre_path: Optional[str]

    ITP_PATH_FORMAT: str = "{molecule_name}_GMX.itp"
    GRO_PATH_FORMAT: str = "{molecule_name}_GMX.gro"
    TOP_PATH_FORMAT: str = "{molecule_name}_GMX.top"
    POSRE_PATH_FORMAT: str = "posre_{molecule_name}.itp"

    @property
    def itp_path(self) -> Optional[str]:
        return f"{self.molecule_name}_GMX.itp" if self.itp else None

    @property
    def gro_path(self) -> Optional[str]:
        return f"{self.molecule_name}_GMX.gro" if self.gro else None

    @property
    def top_path(self) -> Optional[str]:
        return f"{self.molecule_name}_GMX.top" if self.top else None

    @property
    def posre_path(self) -> Optional[str]:
        return f"posre_{self.molecule_name}.itp" if self.posre else None


SOLVENT_ITP_DIR = "preprocessed_output/solvents/solvent_itp"
SOLVENT_MOL2_DIR = "preprocessed_output/solvents/solvent_mol2"
SOLVENT_JSON_PATH = "preprocessed_output/solvents/solvent_atomtypes.json"
POLYMER_STRUCURE_OUTPUT = "preprocessed_output/polymer_structures"

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
