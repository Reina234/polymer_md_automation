from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1", "C=C"])
generator.generate_polymer(5, "rdkit_test3", overwrite=False, save=True)
