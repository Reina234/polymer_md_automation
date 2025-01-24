from mappers.martini_index_generator import MARTINIIndexGenerator
from mappers.martini_map_generator import MARTINIMapGenerator
from rdkit_new.base_polymer_generator import BasePolymerGenerator

from dataclasses import dataclass


class MARTINIMaps:
    """
    Manages the generation of MARTINI .map and .ndx files.
    """

    def __init__(self, polymer: BasePolymerGenerator, file_name: str, output_dir: str):
        self.polymer = polymer
        self.file_name = file_name
        self.output_dir = output_dir

    @property
    def map_file(self) -> str:
        """
        Generates the MARTINI .map file only when accessed.
        """
        return MARTINIMapGenerator(self.polymer).create_map(
            self.file_name, self.output_dir
        )

    @property
    def ndx_file(self) -> str:
        """
        Generates the MARTINI .ndx file only when accessed.
        """
        return MARTINIIndexGenerator(self.polymer).save_index(
            self.file_name, self.output_dir
        )
