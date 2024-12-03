from A_modules.shared.file_processing.file_parsers.pdb_parser import PDBParser
from A_modules.shared.utils import check_file_exists, get_file_contents


pdb_parser = PDBParser()
print(check_file_exists("styrene.pdb"))
print(get_file_contents("styrene.pdb"))
content = PDBParser.read_file("styrene.pdb")
print(pdb_parser.get_atom_coordinates(content))
print(pdb_parser.extract_box_dimensions(content))
