from rdkit_new.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1", "C=C"])
generator.generate_polymer(5, "rdkit_test2", overwrite=False, save=False)
print(generator.votca_bonds)
print(generator.votca_map)
print(generator.votca_angles)
print(generator.sequence)
# print(generator.map)
# print(generator.votca_map)
# print(generator.pycg_map)
mapper = VOTCAMappingGenerator(
    "POLY",
    generator.votca_map,
    generator.votca_bonds,
    generator.votca_angles,
    itp_file_path="alternating_test/POLY_GMX.itp",
)
mapper.save_to_xml("rdkit_test4.xml")
