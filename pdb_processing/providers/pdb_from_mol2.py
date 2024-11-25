import os
import parmed as pmd
from pdb_processing.providers.base_pdb_provider import PDBProvider
from solvent.solvent import Solvent
from typing import Optional
class PDBFromMOL2(PDBProvider):
    """
    Generates a PDB file from a MOL2 file.
    """

    def __init__(self, mol2_file: str, solvent: Solvent, additional_notes: Optional[str] = None):
        default_note = f"Loaded from MOL2 file at: {mol2_file}"
        super().__init__(default_note, solvent, additional_notes)
        self.mol2_file = mol2_file

    def get_pdb(self, output_dir: str, identifier: str = "") -> str:
        """
        Convert the MOL2 file into a PDB file with a consistent naming scheme.
        """
        if not os.path.exists(self.mol2_file):
            raise FileNotFoundError(f"MOL2 file not found: {self.mol2_file}")

        identifier_suffix = f"_{identifier}" if identifier else ""
        output_pdb = os.path.join(output_dir, f"{self._molecule_name}_MOL2{identifier_suffix}.pdb")
        os.makedirs(output_dir, exist_ok=True)

        structure = pmd.load_file(self.mol2_file)
        structure.save(output_pdb, format="PDB")
        print(f"PDB file generated from MOL2: {output_pdb}")
        return output_pdb

    def get_source(self) -> str:
        return "MOL2"