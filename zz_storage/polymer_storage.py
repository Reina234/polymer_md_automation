import os
import json
from pathlib import Path
from typing import Dict, Optional, List
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from collections import defaultdict
from itp_parser import ITPModifier
from config.acpype_config import AcpypeOutputConfig
import re
from data_models.output_types import GromacsPaths
from zz_storage.residue_itp_storage import ResidueITPStorage


class PolymerStorage:
    """Handles storage & retrieval of `n=3` parameterized monomers."""

    def __init__(
        self,
        cache_dir: str = "parameterized_monomers",
        residue_storage: ResidueITPStorage = ResidueITPStorage(),
    ):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.monomer_record_file = self.cache_dir / "monomer_records.json"
        self.monomer_records = self._load_monomer_records()
        self.residue_storage = residue_storage
        self.atomtypes_file = self.cache_dir / "atomtypes.json"
        self.atomtypes = self._load_atomtypes()

    def _load_atomtypes(self):
        """Loads stored atomtypes from disk."""
        if self.atomtypes_file.exists():
            with open(self.atomtypes_file, "r") as f:
                return json.load(f)
        return {}

    def store_atomtypes(self, new_atomtypes):
        """Stores new atomtype values, averaging duplicates."""
        for atom, values in new_atomtypes.items():
            if atom in self.atomtypes:
                # ✅ Average values for duplicates
                self.atomtypes[atom]["mass"] = sum(values["mass"]) / len(values["mass"])
                self.atomtypes[atom]["charge"] = sum(values["charge"]) / len(
                    values["charge"]
                )
                self.atomtypes[atom]["ptype"] = "A"  # Default assumption
                self.atomtypes[atom]["sigma"] = sum(values["sigma"]) / len(
                    values["sigma"]
                )
            else:
                self.atomtypes[atom] = {
                    "mass": sum(values["mass"]) / len(values["mass"]),
                    "charge": sum(values["charge"]) / len(values["charge"]),
                    "ptype": "A",
                    "sigma": sum(values["sigma"]) / len(values["sigma"]),
                }

        # ✅ Save to disk
        with open(self.atomtypes_file, "w") as f:
            json.dump(self.atomtypes, f, indent=4)

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

        if self.monomer_records is None:
            print("⚠️ ERROR: monomer_records is None!")
            return False  # Prevent error

        exists = monomer_smiles in self.monomer_records

        return exists

    def store_monomer_itp(
        self, monomer_smiles: str, residue_smiles_list: List[str], itp_dict: Dict
    ):
        """Stores `.itp` data for each residue SMILES."""
        monomer_file = self._get_monomer_file(monomer_smiles)
        with open(monomer_file, "w") as f:
            json.dump(itp_dict, f, indent=4)

        self.monomer_records[monomer_smiles] = residue_smiles_list
        self._save_monomer_records()

    def load_monomer_itp(self, monomer_smiles: str) -> Optional[Dict]:
        """Loads `.itp` data for a monomer if it exists."""
        if not self.monomer_exists(monomer_smiles):
            return None
        monomer_file = self._get_monomer_file(monomer_smiles)
        with open(monomer_file, "r") as f:
            return json.load(f)

    def parameterize_and_store(
        self,
        monomer_smiles: str,
        polymer_generator: BasePolymerGenerator,
        parametizer: ACPYPEParameterizer = ACPYPEParameterizer(),
    ):
        """Runs `n=3` parameterization & stores `.itp` data if missing."""
        print("!!!!", self.monomer_exists(monomer_smiles))
        if self.monomer_exists(monomer_smiles):
            print("!!!exists")
            return  # Already stored, no need to recompute

        print(f"⚠️ Parameterizing `{monomer_smiles}` (n=3)")

        # **Generate `n=3` polymer**
        generator = type(polymer_generator)([monomer_smiles])
        pdb_path = generator.generate_polymer(num_units=3, output_dir=self.cache_dir)

        # **Run ACPYPE parameterization**
        file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
        gromacs_paths = parametizer.run(
            pdb_path,
            self.cache_dir,
            new_file_name=monomer_smiles,
            acpype_output_config=file_config,
            skip_existing=True,
        )

        # **Extract `.itp` contents**
        itp_sections = self.parse_itp_file(gromacs_paths.itp_path)

        # **Segment `.itp` data based on residues**
        segmented_itp = self._segment_itp_by_residue(generator, itp_sections)

        if not segmented_itp:
            raise ValueError(f"⚠️ ERROR: No residues segmented from {monomer_smiles}!")

        # **Store residue `.itp` separately via `ResidueITPStorage`**
        self.residue_storage.store_residue_itp(segmented_itp)  # ✅ FIXED!

        # **Store `.itp` data at the monomer level**
        residue_smiles_list = list(segmented_itp.keys())
        self.store_monomer_itp(monomer_smiles, residue_smiles_list, segmented_itp)

        # **✅ Update `self.monomer_records` to ensure it's properly recorded**
        self.monomer_records[monomer_smiles] = residue_smiles_list
        self._save_monomer_records()

    def _segment_itp_by_residue(
        self, polymer_generator: BasePolymerGenerator, itp_sections: Dict
    ) -> Dict:
        """Segments `.itp` data based on residue SMILES from `self.sequence`, ensuring correct headers."""

        segmented_itp = defaultdict(
            lambda: {
                "atoms": [],
                "bonds": [],
                "angles": [],
                "dihedrals": [],
                "improper_dihedrals": [],
                "constraints": [],
            }
        )

        # ✅ Define expected headers for each section
        section_headers = {
            "atoms": ["nr", "type", "resnr", "res", "atom", "cgnr", "charge", "mass"],
            "bonds": ["ai", "aj", "funct", "r", "k"],
            "angles": ["ai", "aj", "ak", "funct", "theta"],
            "dihedrals": ["ai", "aj", "ak", "al", "funct", "phase", "kd", "pn"],
            "improper_dihedrals": ["ai", "aj", "ak", "al", "funct"],
            "constraints": ["ai", "aj", "funct"],
        }

        # ✅ Get the ordered sequence of residue SMILES
        residue_sequence = polymer_generator.sequence

        # ✅ Ensure mapping between residue SMILES and their atom indices
        residue_atom_mapping = {
            entry["smiles"]: set(entry["atom_indices"])  # Use a set for fast lookup
            for entry in polymer_generator.votca_map
        }

        # ✅ STEP 1: Assign ATOMS to the correct residue
        for section, content in itp_sections.items():
            for line in content:
                if isinstance(line, str):
                    row_parts = line.split()

                    try:
                        if section == "atoms":
                            row = dict(zip(section_headers["atoms"], row_parts))
                            row["nr"] = int(row["nr"])
                            row["resnr"] = int(row["resnr"])
                            row["cgnr"] = int(row["cgnr"])
                            row["charge"] = float(row["charge"])
                            row["mass"] = float(row["mass"])
                        elif section == "bonds":
                            row = dict(zip(section_headers["bonds"], row_parts))
                            row["ai"], row["aj"], row["funct"] = map(
                                int, [row["ai"], row["aj"], row["funct"]]
                            )
                            row["r"], row["k"] = map(float, [row["r"], row["k"]])
                        elif section == "angles":
                            row = dict(zip(section_headers["angles"], row_parts))
                            row["ai"], row["aj"], row["ak"], row["funct"] = map(
                                int, [row["ai"], row["aj"], row["ak"], row["funct"]]
                            )
                            row["theta"] = float(row["theta"])
                        elif section == "dihedrals":
                            row = dict(zip(section_headers["dihedrals"], row_parts))
                            row["ai"], row["aj"], row["ak"], row["al"], row["funct"] = (
                                map(
                                    int,
                                    [
                                        row["ai"],
                                        row["aj"],
                                        row["ak"],
                                        row["al"],
                                        row["funct"],
                                    ],
                                )
                            )
                            row["phase"], row["kd"], row["pn"] = map(
                                float, [row["phase"], row["kd"], row["pn"]]
                            )
                        elif section == "improper_dihedrals":
                            row = dict(
                                zip(section_headers["improper_dihedrals"], row_parts)
                            )
                            row["ai"], row["aj"], row["ak"], row["al"], row["funct"] = (
                                map(
                                    int,
                                    [
                                        row["ai"],
                                        row["aj"],
                                        row["ak"],
                                        row["al"],
                                        row["funct"],
                                    ],
                                )
                            )
                        elif section == "constraints":
                            row = dict(zip(section_headers["constraints"], row_parts))
                            row["ai"], row["aj"], row["funct"] = map(
                                int, [row["ai"], row["aj"], row["funct"]]
                            )
                        else:
                            continue  # Skip unrecognized sections
                    except (IndexError, ValueError):
                        continue  # Skip invalid lines
                else:
                    row = line  # If it's already a dictionary, use it directly

                # ✅ STEP 2: Assign the row to the correct residue(s)
                atom_indices = {
                    row[key] for key in ["ai", "aj", "ak", "al"] if key in row
                } | {row.get("nr", None)}
                atom_indices.discard(None)  # Remove None values

                matched_residues = set()
                for residue_smiles, indices in residue_atom_mapping.items():
                    if any(idx in indices for idx in atom_indices):
                        matched_residues.add(residue_smiles)

                # ✅ STEP 3: Handle assignments correctly
                if section == "atoms":
                    # ✅ Assign ATOMS to the correct residue
                    for residue_smiles in matched_residues:
                        segmented_itp[residue_smiles]["atoms"].append(row)

                elif section in ["bonds", "angles", "dihedrals", "improper_dihedrals"]:
                    if len(matched_residues) == 1:
                        # ✅ If the entire bond/angle/dihedral is **within one residue**, assign normally
                        residue_smiles = next(iter(matched_residues))
                        segmented_itp[residue_smiles][section].append(row)
                    elif len(matched_residues) > 1:
                        # 🚀 **SPECIAL CASE: Bond/Angle/Dihedral connects multiple residues**
                        residue_list = sorted(
                            matched_residues, key=residue_sequence.index
                        )
                        first_residue = residue_list[0]
                        segmented_itp[first_residue][section].append(
                            row
                        )  # Assign to the first residue in sequence

                else:
                    # ✅ Catch unexpected cases
                    print(f"⚠️ WARNING: Skipping unexpected entry in {section}: {row}")

        return segmented_itp

    @staticmethod
    def parse_itp_file(itp_path):
        """
        Parses a GROMACS `.itp` file and extracts all sections, distinguishing proper and improper dihedrals.
        """
        print(itp_path)
        sections = defaultdict(list)
        current_section = None
        is_improper = False

        with open(itp_path, "r") as f:
            for line in f:
                line = line.strip()

                if not line or line.startswith(";"):
                    continue

                section_match = re.match(r"\[(.+)\]", line)
                if section_match:
                    current_section = section_match.group(1).strip()
                    is_improper = False  # Reset improper flag

                    # If it's dihedrals, check if it's improper
                    if current_section == "dihedrals":
                        line_after = f.readline().strip()
                        if "improper" in line_after.lower():
                            current_section = "improper_dihedrals"
                        else:
                            sections[current_section].append(
                                line_after
                            )  # Retain first data row

                    continue

                if current_section:
                    sections[current_section].append(line)

        return sections
