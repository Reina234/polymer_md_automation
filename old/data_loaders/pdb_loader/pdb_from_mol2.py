# providers/pdb_mol2_provider.py
import parmed as pmd
from .base_pdb_provider import PDBProvider


class PDBFromMOL2(PDBProvider):
    def __init__(self, mol2_file: str):
        self.mol2_file = mol2_file

    def get_pdb(self, output_file: str) -> None:
        """
        Convert a MOL2 file to a PDB file.
        """
        structure = pmd.load_file(self.mol2_file)
        structure.save(output_file, format="PDB")
        print(f"PDB file converted from MOL2 and saved to: {output_file}")
