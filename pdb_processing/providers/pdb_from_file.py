import shutil
from typing import Optional
from pdb_processing.providers.base_pdb_provider import PDBProvider
import os
from solvent.solvent import Solvent
class PDBFromFile(PDBProvider):
    """
    Loads a PDB file from an existing source.
    """

    def __init__(self, pdb_file: str, solvent: Solvent, additional_notes: Optional[str] = None):
        default_note = f"Loaded from existing PDB file at: {pdb_file}"
        super().__init__(default_note, solvent, additional_notes)
        self.pdb_file = pdb_file

    def get_pdb(self, output_dir: str, identifier: str = "") -> str:
        """
        Copy the PDB file to the output directory with a consistent naming scheme.
        """

        if not os.path.exists(self.pdb_file):
            raise FileNotFoundError(f"PDB file not found: {self.pdb_file}")

        identifier_suffix = f"_{identifier}" if identifier else ""
        output_pdb = os.path.join(output_dir, f"{self._molecule_name}_File{identifier_suffix}.pdb")
        os.makedirs(output_dir, exist_ok=True)

        shutil.copy(self.pdb_file, output_pdb)
        print(f"PDB file loaded from existing file: {output_pdb}")
        return output_pdb

    def get_source(self) -> str:
        return "Downloaded File"
