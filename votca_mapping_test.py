from rdkit_new.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 3, "rdkit_test2", overwrite=False, save=False)
print(generator.votca_bonds)
print(generator.votca_map)
print(generator.votca_angles)
# print(generator.map)
# print(generator.votca_map)
# print(generator.pycg_map)
mapper = VOTCAMappingGenerator(
    "POLY",
    generator.votca_map,
    generator.votca_bonds,
    generator.votca_angles,
    itp_file_path="1_22_test/c=cc1ccccc1_3/POLY_GMX.itp",
)
mapper.save_to_xml("rdkit_test3.xml")
