import subprocess
from typing import Optional
from pdb_processing.providers.base_pdb_provider import PDBProvider
from pdb_processing.utils import validate_pdb_lightweight


class PDBManager:
    """
    Manages PDB generation/loading and metadata.
    """

    def __init__(self, provider: PDBProvider):
        self.provider = provider
        self._pdb_path: Optional[str] = None
        self.solvent = provider.solvent
        #NOTE: is self.solvent needed?

    def prepare_pdb(self, output_dir: str, identifier: str = "") -> str:
        """
        Generate or load a PDB file using the assigned provider and validate it.

        Args:
            output_dir (str): Directory where the PDB file will be saved.
            identifier (str, optional): An identifier for the PDB file. Defaults to "".

        Returns:
            str: The validated PDB file path.

        Raises:
            ValueError: If the PDB validation fails.
        """
        # Generate or load the PDB file
        self._pdb_path = self.provider.get_pdb(output_dir, identifier)

        # Validate the PDB file
        if self._validate_pdb(self._pdb_path):
            raise ValueError(f"PDB validation failed for file: {self._pdb_path}. Please update the PDB file.")

        print("Validation successful. PDB file is ready.")
        return self._pdb_path

    def _validate_pdb(self, pdb_path: str) -> bool:
        """
        Validates the PDB file using pdb-tools. Returns True if validation fails.
        """
        try:
            validate_pdb_lightweight(pdb_path)
            print(f"[+] Validation Successful: {pdb_path}")
            return False
        except ValueError as e:
            print(f"[!] Validation Failed for {pdb_path}:\n{e}")
            return True

#NOTE: is this needed? 
    @property
    def metadata(self) -> str:
        """
        Retrieve metadata about the PDB file, including notes.
        """
        return self.provider.metadata

    def update_additional_notes(self, notes: str) -> None:
        """
        Update the additional notes for the PDB provider.
        """
        print("Current notes:", self.provider.additional_notes)
    
        # Prompt user to confirm update
        overwrite = input("Do you want to update the notes? [y/n]: ").strip().lower()
        if overwrite == "no":
            print("[!] No changes made to additional notes. Returning.")
            return
        elif overwrite != "yes":
            raise ValueError("Invalid input. Please enter 'y' for yes or 'n' for no.")
    
        # Perform the update
        self.provider.additional_notes = notes
        print(f"[+] Additional notes updated to: {notes}")

    @property
    def pdb_path(self) -> str:
        """
        Retrieve the PDB path with validation.
        """
        if not self._pdb_path:
            raise ValueError("PDB file has not been prepared yet. Please call .prepare_pdb() first.")
        return self._pdb_path