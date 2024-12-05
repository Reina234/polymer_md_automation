from A_modules.shared.file_processing.file_parsers.pdb_parser import PDBParser
from A_modules.shared.utils import check_file_exists, get_file_contents


pdb_parser = PDBParser()
check_file_exists("styrene.pdb")
get_file_contents("styrene.pdb")
content = PDBParser.read_file("styrene.pdb")

print(PDBParser.get_atom_coordinates(PDBParser.extract_atoms(content)))
