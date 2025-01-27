from modules.rdkit.polymer_builders.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from rdkit import Chem
from modules.cg_mappers.martini_index_generator import MARTINIIndexGenerator

from modules.cg_mappers.votca_map_generator import VOTCAMapGenerator
from modules.cg_mappers.martini_map_generator import MARTINIMapGenerator
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from modules.cg_mappers.pycgtool_map_generator import PyCGToolMapGenerator
from modules.cg_mappers.open_mscg_map_generator import OpenMSCGMapGenerator

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1", "C=C"])
generator.generate_polymer(10, "rdkit_test2", overwrite=False, save=False)
print(generator.cg_map)
