from rdkit_new.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator
from rdkit import Chem
from mappers.martini_index_generator import MARTINIIndexGenerator

from mappers.votca_map_generator import VOTCAMapGenerator
from mappers.martini_map_generator import MARTINIMapGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from mappers.pycgtool_map_generator import PyCGToolMapGenerator
from mappers.open_mscg_map_generator import OpenMSCGMapGenerator

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1", "C=C"])
generator.generate_polymer(10, "rdkit_test2", overwrite=False, save=False)
print(generator.cg_map)
