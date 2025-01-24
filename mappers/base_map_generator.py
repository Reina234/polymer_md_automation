import xml.etree.ElementTree as ET
from typing import List, Dict
import xml.etree.ElementTree as ET
from typing import List, Tuple, Dict, Optional
import xml.dom.minidom
import pandas as pd
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from modules.shared.utils.file_utils import check_directory_exists
from abc import ABC, abstractmethod
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaseMapGenerator(ABC):
    """
    Generates a VOTCA-compatible XML mapping file for coarse-grained simulations.
    """

    map_extension: str

    def __init__(
        self,
        polymer: BasePolymerGenerator,
    ):
        self.bead_mappings = polymer.cg_map
        self.bonds = polymer.cg_bonds if polymer.cg_bonds else []
        self.angles = polymer.cg_angles if polymer.cg_angles else []

    def _generate_filepath(self, filename: str, output_dir) -> str:
        check_directory_exists(output_dir, make_dirs=True)
        return f"{output_dir}/{filename}.{self.map_extension}"

    @abstractmethod
    def _generate_mapping(self) -> str:
        pass

    def _save_mapping(
        self, content, filename: str, output_dir: Optional[str] = None
    ) -> str:
        file_path = self._generate_filepath(filename, output_dir)
        with open(file_path, "w") as f:
            f.write(content)
        return file_path

    def create_map(self, filename: str, output_dir: Optional[str] = None) -> str:
        content = self._generate_mapping()
        return self._save_mapping(content, filename, output_dir)

    @staticmethod
    def reformat_atom_0(atom_name: str) -> str:
        """
        Removes the '0' from the first occurrence of an atom type (e.g., 'C0' -> 'C', 'H0' -> 'H').
        Ensures that 'C10', 'H10', etc., remain unchanged.

        Args:
            atom_name (str): The raw atom name (e.g., 'C0', 'H0', 'C10').

        Returns:
            str: The formatted atom name (e.g., 'C' instead of 'C0').
        """
        match = re.match(r"^([A-Z]+)(0)$", atom_name)  # Match only if it's AtomType + 0
        if match:
            return match.group(1)  # Return only the atom type (e.g., 'C' for 'C0')

        return atom_name
