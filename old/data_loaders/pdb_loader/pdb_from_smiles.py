# providers/pdb_rdkit_provider.py
from openff.toolkit.topology import Molecule
from .base_pdb_provider import PDBProvider


class PDBFromSmiles(PDBProvider):
    def __init__(self, smiles: str):
        self.smiles = smiles

    def get_pdb(self, output_file: str) -> None:
        """
        Generate a PDB file from a SMILES string using RDKit.
        """
        molecule = Molecule.from_smiles(self.smiles)
        molecule.generate_conformers()
        molecule.to_file(output_file, file_format="pdb")
        print(f"PDB file generated from SMILES and saved to: {output_file}")
