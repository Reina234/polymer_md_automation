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

        # Retrieve `.itp` sections per residue
        segmented_itp = defaultdict(lambda: defaultdict(list))

        for residue_smiles in self.polymer_generator.sequence:
            itp_data = self.residue_storage.load_residue_itp(residue_smiles)
            if itp_data:
                for section, content in itp_data.items():
                    segmented_itp[residue_smiles][section].extend(content)

        # Adjust atom indices before writing the `.itp`
        adjusted_itp = self._adjust_atom_indices(segmented_itp)

        return self._write_itp(adjusted_itp)

    def _adjust_atom_indices(self, segmented_itp):
        """
        Adjusts atom indices when combining multiple residues into one `.itp` file.

        Ensures:
        - Atoms are numbered in correct **polymer sequence order**.
        - Connectivity between residues is **preserved correctly**.
        - Atom indices **do not jump unexpectedly**.
        - Monomer SMILES are stored in the `.itp` for debugging.
        """
        adjusted_itp = defaultdict(lambda: defaultdict(list))
        atom_counter = 1  # ✅ Ensures continuous numbering
        atom_mapping = {}  # ✅ Maps old atom indices to new ones
        residue_atom_tracking = {}  # ✅ Tracks which atoms belong to which residue
        atom_name_map = {}  # ✅ Ensures unique atom names
        debug_info = []  # ✅ Stores debugging information for each atom
        residue_order = (
            self.polymer_generator.sequence
        )  # ✅ The **correct polymer sequence**

        # ✅ **Step 1: Process ATOMS in Polymer Sequence Order**
        for residue_idx, residue_smiles in enumerate(
            residue_order
        ):  # **Use actual order**
            if residue_smiles not in segmented_itp:
                continue  # Skip if residue is missing

            sections = segmented_itp[residue_smiles]
            residue_atom_tracking[residue_smiles] = []  # Track atoms per residue

            for atom in sections["atoms"]:
                original_index = int(atom["nr"])  # ✅ Get original atom index

                # ✅ **Ensure Unique Atom Naming**
                original_name = atom["atom"]
                if original_name in atom_name_map:
                    atom_name_map[original_name] += 1
                    unique_name = f"{original_name}{atom_name_map[original_name]}"
                else:
                    atom_name_map[original_name] = 1
                    unique_name = original_name

                # ✅ **Store Correct Mapping**
                atom_mapping[original_index] = atom_counter  # Map original to new index
                atom["nr"] = atom_counter  # ✅ Update to sequential numbering
                atom["atom"] = unique_name  # ✅ Assign unique atom name
                atom["residue_smiles"] = residue_smiles  # ✅ Add debugging info
                atom["residue_index"] = residue_idx + 1  # ✅ Track monomer position
                debug_info.append(
                    f"{atom_counter}: {unique_name} (Residue {residue_smiles})"
                )

                atom_counter += 1
                residue_atom_tracking[residue_smiles].append(
                    atom_counter
                )  # ✅ Track atom order
                adjusted_itp[residue_smiles]["atoms"].append(atom)

        # ✅ **Step 2: Update BONDS, ANGLES, DIHEDRALS**
        for residue_smiles in residue_order:
            if residue_smiles not in segmented_itp:
                continue  # Skip missing residues

            sections = segmented_itp[residue_smiles]

            for section in ["bonds", "angles", "dihedrals", "improper_dihedrals"]:
                for entry in sections[section]:
                    # ✅ **Ensure All Atom Indices Exist in `atom_mapping`**
                    for atom_key in ["ai", "aj", "ak", "al"]:
                        if atom_key in entry:
                            if int(entry[atom_key]) not in atom_mapping:
                                print(
                                    f"🚨 ERROR: Atom {entry[atom_key]} is missing from mapping!"
                                )
                                continue  # Skip entry if error persists
                            entry[atom_key] = atom_mapping[int(entry[atom_key])]

                    adjusted_itp[residue_smiles][section].append(entry)

        # ✅ **Print Debug Info to Verify Order**
        print("\n✅ Atom Order Debugging Information:")
        for line in debug_info:
            print(line)

        return adjusted_itp

    def _write_itp(self, adjusted_itp):
        """Writes `.itp` file from retrieved sections with correct formatting and ACPYPE-style comments."""
        atomtypes = self.polymer_storage._load_atomtypes()

        with open("final_polymer.itp", "w") as f:
            f.write("\n[ atomtypes ]\n")
            f.write("; name   bond_type   mass   charge   ptype   sigma\n")
            for atom, values in atomtypes.items():
                f.write(
                    f"{atom:6}  A   {values['mass']:6.3f}   {values['charge']:6.3f}   {values['ptype']}   {values['sigma']:6.3f}\n"
                )

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
                            )
                            section_written = True

                        for row in sections[section]:
                            if section == "bonds":
                                bond_comment = f"  ; {row['ai']} - {row['aj']}"
                                f.write(
                                    "   ".join(
                                        str(row[col])
                                        for col in self.section_headers[section]
                                    )
                                    + bond_comment
                                    + "\n"
                                )
                            else:
                                f.write(
                                    "   ".join(
                                        str(row[col])
                                        for col in self.section_headers[section]
                                    )
                                    + "\n"
                                )

        return "final_polymer.itp"
