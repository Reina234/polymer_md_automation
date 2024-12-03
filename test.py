from A_modules.shared.file_processing.file_parsers.pdb_parser import PDBParser
from A_modules.shared.utils import check_file_exists


pdb_parser = PDBParser()
content = PDBParser.read_file("styrene.pdb")
print(PDBParser().extract_atoms(content))
