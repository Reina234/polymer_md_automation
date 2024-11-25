from rdkit import Chem
from rdkit.Chem import AllChem
import os
from typing import Optional
from pdb_processing.providers.base_pdb_provider import PDBProvider
from solvent.solvent import Solvent
from pdb_processing.utils import handle_existing_file

class PDBFromSmilesRDKit(PDBProvider):
    """
    Generates a PDB file from a SMILES string using RDKit.
    """

    source = "RDkit"

    def __init__(self, smiles: str, solvent: Solvent, additional_notes: Optional[str] = None):
        default_note = f"Generated from SMILES: {smiles} using RDkit"
        super().__init__(default_note, solvent, additional_notes)
        self.smiles = smiles

    def get_pdb(self, output_dir: str, identifier: str = "") -> str:
        """
        Generate a PDB file from a SMILES string with a consistent naming scheme.
        """
        molecule = Chem.MolFromSmiles(self.smiles)
        if molecule is None:
            raise ValueError(f"Invalid SMILES string: {self.smiles}")

        molecule = Chem.AddHs(molecule)
        AllChem.EmbedMolecule(molecule)
        AllChem.UFFOptimizeMolecule(molecule)

        identifier_suffix = f"_{identifier}" if identifier else ""
        output_pdb = os.path.join(output_dir, f"{self._molecule_name}_{self.source}{identifier_suffix}.pdb")


        os.makedirs(output_dir, exist_ok=True)

        existing_file = handle_existing_file(output_pdb)
        if existing_file:
            return existing_file

        Chem.MolToPDBFile(molecule, output_pdb)
        print(f"[+] PDB file generated from RDKit: {output_pdb}")    

        return output_pdb


