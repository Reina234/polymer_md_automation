from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from rdkit import Chem
from zz_storage.gromacs_file_builder import GromacsFileBuilder
from modules.rdkit.votca_mapping_generator import VOTCAMappingGenerator
from modules.gromacs.parsers.gromacs_parser import GromacsParser
from modules.gromacs.parsers.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple
from zz_storage.residue_itp_storage import itp_storage
from zz_storage.polymer_storage import PolymerStorage

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1", "C=C"])
generator.generate_polymer(5, "rdkit_test3", overwrite=True, save=True)

gromacs = GromacsFileBuilder(generator, PolymerStorage(), itp_storage)
gromacs.build_gromacs_files()
