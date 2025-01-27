from modules.rdkit.polymer_builders.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem

from modules.rdkit.votca_mapping_generator import VOTCAMappingGenerator
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 3, "rdkit_test2", overwrite=False, save=False)
print(generator.cg_bonds)
print(generator.cg_map)
print(generator.cg_angles)
# print(generator.map)
# print(generator.votca_map)
# print(generator.pycg_map)
mapper = VOTCAMappingGenerator(
    "POLY",
    generator.cg_map,
    generator.cg_bonds,
    generator.cg_angles,
    itp_file_path="1_22_test/c=cc1ccccc1_3/POLY_GMX.itp",
)
mapper.save_to_xml("rdkit_test3.xml")
