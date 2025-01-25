from typing import List
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.utils.file_utils import calculate_minimum_box_size_from_df
from modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler


def add_box_dim(
    gro_file, padding: float = 0.1, parser=GromacsParser(), gro_handler=GroHandler()
) -> List[float]:
    sections = parser.parse(gro_file)
    first_key = next(iter(sections))  # Get the first key
    gro_section = sections[first_key]
    gro_handler.process(gro_section)
    box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding)
    gro_handler.box_dimensions = box_size
    sections[first_key] = gro_handler.export()

    parser.export(sections, gro_file)

    return box_size
