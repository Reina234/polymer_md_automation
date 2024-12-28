import os
from enum import Enum
import dataclasses
from typing import Optional

TEMP_DIR = "temp"
LOG_DIR = "logs"
TOPOL_NAME = "topol.top"
CACHE_DIR = "cache"
MDP_CACHE_DIR = os.path.join(CACHE_DIR, "mdp_cache")

# SOLVENT_ITP_DIR = "preprocessed_output/solvents/solvent_itp"
# SOLVENT_MOL2_DIR = "preprocessed_output/solvents/solvent_mol2"
# SOLVENT_JSON_PATH = "preprocessed_output/solvents/solvent_atomtypes.json"


# ACPYPE_SOLVENT_OUTPUT_SUBDIR = "acpype_solvent_output"

ACPYPE_POLYMER_NAME = "POLY"
