from modules.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
import os
from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from rdkit import Chem

from modules.rdkit.votca_mapping_generator import VOTCAMappingGenerator
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator.generate_polymer(10, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

path = os.path.abspath("rdkit_test3/c=cc1ccccc1_10.pdb")
EditconfPDBtoGROConverter().run(path, "rdkit_test3")
