import os
import json
from pathlib import Path
from typing import Dict, Optional, List
from rdkit import Chem


class PolymerStorage:
    """Handles storage & retrieval of `n=3` parameterized monomers."""

    def __init__(self, cache_dir: str = "parameterized_monomers"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.monomer_record_file = self.cache_dir / "monomer_records.json"
        self.monomer_records = self._load_monomer_records()

    def _get_monomer_file(self, monomer_smiles: str) -> Path:
        """Returns the storage path for a monomer's `.itp` data."""
        return self.cache_dir / f"{monomer_smiles}.json"

    def _load_monomer_records(self) -> Dict[str, List[str]]:
        """Loads the record of processed monomers & their residue SMILES."""
        if self.monomer_record_file.exists():
            with open(self.monomer_record_file, "r") as f:
                return json.load(f)
        return {}

    def _save_monomer_records(self):
        """Saves the monomer records to disk."""
        with open(self.monomer_record_file, "w") as f:
            json.dump(self.monomer_records, f, indent=4)

    def monomer_exists(self, monomer_smiles: str) -> bool:
        """Checks if `n=3` parameterization exists for a given monomer."""
        return monomer_smiles in self.monomer_records

    def store_monomer_itp(
        self, monomer_smiles: str, residue_smiles_list: List[str], itp_dict: Dict
    ):
        """Stores `.itp` data for each residue SMILES."""
        monomer_file = self._get_monomer_file(monomer_smiles)
        with open(monomer_file, "w") as f:
            json.dump(itp_dict, f, indent=4)

        # Store residue SMILES reference for this monomer
        self.monomer_records[monomer_smiles] = residue_smiles_list
        self._save_monomer_records()

    def load_monomer_itp(self, monomer_smiles: str) -> Optional[Dict]:
        """Loads `.itp` data for a monomer if it exists."""
        if not self.monomer_exists(monomer_smiles):
            00000000000000000000000000000
            return None
        monomer_file = self._get_monomer_file(monomer_smiles)
        with open(monomer_file, "r") as f:
            return json.load(f)

    def remove_monomer_itp(self, monomer_smiles: str):
        """Removes stored `.itp` data for a monomer & its associated residue SMILES."""
        if monomer_smiles not in self.monomer_records:
            return

        # Remove `.itp` file
        monomer_file = self._get_monomer_file(monomer_smiles)
        if monomer_file.exists():
            monomer_file.unlink()

        # Remove monomer record
        del self.monomer_records[monomer_smiles]
        self._save_monomer_records()
