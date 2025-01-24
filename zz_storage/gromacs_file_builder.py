import json
import os
from collections import defaultdict
from typing import Dict, List
from pathlib import Path
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from zz_storage.polymer_storage import PolymerStorage
from zz_storage.residue_itp_storage import ResidueITPStorage


class GromacsFileBuilder:
    """Builds `.itp`, `.gro`, `.top` for a **long polymer**."""

    def __init__(
        self,
        polymer_generator: BasePolymerGenerator,
        polymer_storage: PolymerStorage,
        residue_storage: ResidueITPStorage,
    ):
        self.polymer_generator = polymer_generator
        self.polymer_storage = polymer_storage
        self.residue_storage = residue_storage
        self.section_headers = {
            "atoms": ["nr", "type", "resnr", "res", "atom", "cgnr", "charge", "mass"],
            "bonds": ["ai", "aj", "funct", "r", "k"],
            "angles": ["ai", "aj", "ak", "funct", "theta"],
            "dihedrals": ["ai", "aj", "ak", "al", "funct", "phase", "kd", "pn"],
            "improper_dihedrals": ["ai", "aj", "ak", "al", "funct"],
            "constraints": ["ai", "aj", "funct"],
        }

    def build_gromacs_files(self):
        """Ensures all monomers exist, retrieves `.itp` sections per **residue SMILES**, and builds final `.itp`."""
        for monomer_smiles in self.polymer_generator.monomer_smiles_list:
            if not self.polymer_storage.monomer_exists(monomer_smiles):
                self.polymer_storage.parameterize_and_store(
                    monomer_smiles, self.polymer_generator
                )
        print("!!!!!!!!!!!!!")
        # **Retrieve `.itp` sections per residue**
        segmented_itp = defaultdict(lambda: defaultdict(list))

        for residue_smiles in self.polymer_generator.sequence:
            itp_data = self.residue_storage.load_residue_itp(residue_smiles)
            if itp_data:
                for section, content in itp_data.items():
                    segmented_itp[residue_smiles][section].extend(content)

        # **Adjust atom indices properly before writing the `.itp`**
        adjusted_itp = self._adjust_atom_indices(segmented_itp)

        return self._write_itp(adjusted_itp)

    def _adjust_atom_indices(self, segmented_itp):
        """
        Adjusts atom indices when combining multiple residues into one `.itp` file.

        - Ensures atom numbers in bonds, angles, dihedrals are updated sequentially.
        - Assigns unique atom names while preserving original identity.
        """
        adjusted_itp = defaultdict(lambda: defaultdict(list))
        atom_counter = 1
        atom_mapping = {}
        atom_name_map = {}  # 🚀 Track unique atom names to avoid duplicates

        # ✅ STEP 1: Update ATOMS section & store atom mapping
        for residue_smiles, sections in segmented_itp.items():
            for atom in sections["atoms"]:
                original_index = int(atom["nr"])  # ✅ Fix: Access dictionary key

                # 🚀 Ensure unique atom names (avoid duplicates)
                original_name = atom["atom"]
                if original_name in atom_name_map:
                    atom_name_map[original_name] += 1
                    unique_name = f"{original_name}{atom_name_map[original_name]}"
                else:
                    atom_name_map[original_name] = 1
                    unique_name = original_name

                # 🚀 Store correct mapping
                atom_mapping[original_index] = atom_counter
                atom["nr"] = atom_counter  # ✅ Update atom index correctly
                atom["atom"] = unique_name  # ✅ Assign unique atom name
                atom_counter += 1

                adjusted_itp[residue_smiles]["atoms"].append(atom)

        # ✅ STEP 2: Update BONDS, ANGLES, DIHEDRALS, IMPROPER DIHEDRALS
        for residue_smiles, sections in segmented_itp.items():
            for section in ["bonds", "angles", "dihedrals", "improper_dihedrals"]:
                for entry in sections[section]:
                    # ✅ Ensure all atom indices exist in `atom_mapping`
                    for atom_key in ["ai", "aj", "ak", "al"]:
                        if atom_key in entry:
                            if int(entry[atom_key]) not in atom_mapping:
                                print(f"🚨 ERROR: Atom {entry[atom_key]} is missing!")
                                continue  # 🚨 Skip entry if the issue persists
                            entry[atom_key] = atom_mapping[int(entry[atom_key])]

                    adjusted_itp[residue_smiles][section].append(entry)

        return adjusted_itp

    def _write_itp(self, adjusted_itp):
        """Writes `.itp` file from retrieved sections with correct formatting and headers."""

        # ✅ Load atomtypes before writing
        atomtypes = self.polymer_storage._load_atomtypes()

        with open("final_polymer.itp", "w") as f:
            # ✅ **Write ATOMTYPES first**
            f.write("\n[ atomtypes ]\n")
            f.write("; name   bond_type   mass   charge   ptype   sigma\n")
            for atom, values in atomtypes.items():
                f.write(
                    f"{atom:6}  A   {values['mass']:6.3f}   {values['charge']:6.3f}   {values['ptype']}   {values['sigma']:6.3f}\n"
                )

            # ✅ **Write all sections**
            for section in [
                "atoms",
                "bonds",
                "angles",
                "dihedrals",
                "improper_dihedrals",
                "constraints",
            ]:
                section_written = False
                for residue_smiles, sections in adjusted_itp.items():
                    if section in sections:
                        if not section_written:
                            f.write(f"\n[{section}]\n")
                            f.write(
                                "; " + "   ".join(self.section_headers[section]) + "\n"
                            )  # ✅ Header row
                            section_written = True

                        for row in sections[section]:
                            if isinstance(row, dict):
                                f.write(
                                    "   ".join(
                                        str(row[col])
                                        for col in self.section_headers[section]
                                    )
                                    + "\n"
                                )
                            else:
                                f.write(str(row) + "\n")

        return "final_polymer.itp"


def _write_itp(self, adjusted_itp):
    """Writes `.itp` file from retrieved sections, ensuring proper formatting."""
    # ✅ Load atomtypes before writing
    atomtypes = self.polymer_storage._load_atomtypes()

    with open("final_polymer.itp", "w") as f:
        # ✅ **Write ATOMTYPES first**
        f.write("\n[ atomtypes ]\n")
        f.write("; name   bond_type   mass   charge   ptype   sigma\n")
        for atom, values in atomtypes.items():
            f.write(
                f"{atom:6}  A   {values['mass']:6.3f}   {values['charge']:6.3f}   {values['ptype']}   {values['sigma']:6.3f}\n"
            )

        for section in ["atoms", "bonds", "angles", "dihedrals", "improper_dihedrals"]:
            if section in adjusted_itp:
                f.write(f"\n[{section}]\n")

                # ✅ Add Commented Headers for Readability
                if section == "atoms":
                    f.write("; nr   type   resnr   res   atom   cgnr   charge   mass\n")
                elif section == "bonds":
                    f.write("; ai   aj   funct   r   k\n")
                elif section == "angles":
                    f.write("; ai   aj   ak   funct   theta\n")
                elif section == "dihedrals":
                    f.write("; ai   aj   ak   al   func   phase   kd   pn\n")

                # ✅ Write the Section Content
                for row in adjusted_itp[section]:
                    f.write(" ".join(map(str, row.values())) + "\n")

    return "final_polymer.itp"
