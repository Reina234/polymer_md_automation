import os
from enum import Enum
import dataclasses
from typing import Optional

TEMP_DIR = "temp"
LOG_DIR = "logs"
TOPOL_NAME = "topol.top"
CACHE_DIR = "cache"
PREPROCESSED_DIR = "preprocessed"
PARAMETERISED_POLYMER_DIR = os.path.join(PREPROCESSED_DIR, "parameterised_polymers")
EQUILIBRIATED_SOLVENT_BOX_DIR = os.path.join(
    PREPROCESSED_DIR, "equilibriated_solvent_boxes"
)
MDP_CACHE_DIR = os.path.join(CACHE_DIR, "mdp_cache")

# SOLVENT_ITP_DIR = "preprocessed_output/solvents/solvent_itp"
# SOLVENT_MOL2_DIR = "preprocessed_output/solvents/solvent_mol2"
# SOLVENT_JSON_PATH = "preprocessed_output/solvents/solvent_atomtypes.json"


# ACPYPE_SOLVENT_OUTPUT_SUBDIR = "acpype_solvent_output"

ACPYPE_POLYMER_NAME = "POLY"
PACKMOL_TEMPLATE_DIR = "modules/shared/packmol/templates"

EQUILIBRIATED_OUTPUTS_SUBDIR = "equilibriated_outputs"
