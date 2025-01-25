from rdkit_new.homopolymer_generator import HomopolymerGenerator
from config.paths import TEMP_DIR
from config.acpype_config import AcpypeOutputConfig
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from modules.shared.file_conversion.converter_factory import ConverterFactory
from modules.shared.utils.file_utils import copy_file, delete_directory
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from data_models.output_types import GromacsPaths
from config.paths import PARAMETERISED_POLYMER_DIR
import os

# NOTE: I think I can use just a copolymer generator, and create a homopolymer by passing in two of the same monomer units


def parameterize(
    polymer_pdb: str,
    output_dir: str = PARAMETERISED_POLYMER_DIR,
    temp_dir: str = TEMP_DIR,
    keep_pdb: bool = False,
    keep_mol2: bool = False,
    polymer_name: str = "POLY",
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    converter_factory: ConverterFactory = ConverterFactory(),
    cleanup: bool = True,
) -> GromacsPaths:

    if keep_pdb:
        polymer_pdb = copy_file(polymer_pdb, output_dir)

    if keep_mol2:
        mol2_output_dir = output_dir
    else:
        mol2_output_dir = temp_dir
    mol2_converter = converter_factory.get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(polymer_pdb, mol2_output_dir)
    parameterizer = parameterizer(acpype_molecule_name=polymer_name)
    file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
    parametized_polymer = parameterizer.run(mol2_file, output_dir, file_config)
    if cleanup:
        delete_directory(temp_dir, confirm=False)
    return parametized_polymer


from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator.generate_polymer(10, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

import numpy as np
import MDAnalysis as mda
from collections import defaultdict
