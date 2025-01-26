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

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
generator.generate_polymer(3, "rdkit_test2", overwrite=False, save=False)
from mappers.open_mscg_map_generator import OpenMSCGMapGenerator

open_msg_map_generator = OpenMSCGMapGenerator(generator, "temp/production.gro")
open_msg_map_generator.create_map("test_open_msg", "zzz")

print(generator.cg_bonds)
print(generator.cg_map)
print(generator.cg_angles)
print(generator.sequence)
# print(generator.map)
# print(generator.votca_map)
# print(generator.pycg_map)
mapper = VOTCAMapGenerator(
    generator,
    itp_file_path="alternating_test/POLY_GMX.itp",
)
mapper.molecule_name = "UNL"
mapper.create_map("test", "zzz")

mapper_martini = MARTINIMapGenerator(generator)
mapper_martini.create_map("test_martini", "zzz")

mapper_index = MARTINIIndexGenerator(generator)
mapper_index.create_map("test_index", "zzz")
from data_models.output_types import MARTINIMaps

martini = MARTINIMaps(generator, "test_all", "zzz")
print(martini.map_file, martini.ndx_file)

pycgtool = PyCGToolMapGenerator(generator)
pycgtool.add_solvent_to_map("1_5_test/hexane/solvent.itp")
pycgtool.create_map("test_pycgtool", "zzz")

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
generator.generate_polymer(10, "rdkit_test2", overwrite=False, save=False)
pycgtool = PyCGToolMapGenerator(generator)
pycgtool.add_solvent_to_map("1_5_test/hexane/solvent.itp")
pycgtool.create_map("test_pycgtool", "zzz", start_index=0)
