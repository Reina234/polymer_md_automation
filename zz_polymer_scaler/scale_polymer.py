import pandas as pd
import os
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import subprocess


class ExpandedPolymerGenerator:
    def __init__(self, base_polymer, num_units, itp_path):
        """
        Expands an existing polymer while keeping the alternating sequence intact.

        Args:
        - base_polymer: Instance of BasePolymerGenerator (from RDKit).
        - num_units (int): Number of monomer units to add.
        - itp_path (str): Path to the .itp file to update.
        """
        self.base_polymer = base_polymer
        self.num_units = num_units
        self.itp_path = itp_path
        self.monomer_smiles = self.base_polymer.monomer_smiles  # List of monomer SMILES
        self.sequence = self._extract_sequence()

    def _extract_sequence(self):
        """
        Extracts the alternating sequence from the n=5 polymer.
        Assumes monomers in `self.monomer_smiles` appear in the correct order.
        """
        return self.monomer_smiles[:-1]  # Exclude terminal monomer

    def _get_next_monomers(self):
        """
        Determines which monomers to add based on the sequence pattern.
        """
        repeat_block = self.sequence  # Middle section of the polymer
        repeat_length = len(repeat_block)

        new_monomers = []
        for i in range(self.num_units):
            next_monomer = repeat_block[i % repeat_length]  # Keep pattern consistent
            new_monomers.append(next_monomer)

        return new_monomers

    def expand_polymer(self, output_itp, output_pdb, output_gro):
        """
        Expands the polymer and writes the new topology.

        Args:
        - output_itp (str): Output .itp file.
        - output_pdb (str): Output .pdb file.
        - output_gro (str): Output .gro file.
        """
        new_monomers = self._get_next_monomers()
        print(f"Adding monomers: {new_monomers}")

        # Load existing topology
        existing_df = convert_itp_to_df(self.itp_path, section="atoms")
        atom_offset = existing_df.shape[0]

        # --- Step 1: Generate New Polymer ---
        polymer_mol = Chem.MolFromSmiles(".".join(self.sequence + new_monomers))
        AllChem.EmbedMolecule(polymer_mol)

        # --- Step 2: Generate New .itp File ---
        with open(output_itp, "w") as f:
            f.write("[ atoms ]\n")
            for idx, monomer in enumerate(new_monomers, start=atom_offset):
                f.write(f"{idx} {monomer} {idx+1}\n")

        print(f"Saved updated topology to {output_itp}")

        # --- Step 3: Save as .pdb ---
        Chem.MolToPDBFile(polymer_mol, output_pdb)
        print(f"Saved expanded polymer structure to {output_pdb}")

        # --- Step 4: Convert to .gro ---
        subprocess.run(
            ["gmx", "editconf", "-f", output_pdb, "-o", output_gro], check=True
        )
        print(f"Converted {output_pdb} to {output_gro}")


# --- Example Usage ---
polymer_generator = BasePolymerGenerator()  # Load existing polymer class
expander = ExpandedPolymerGenerator(
    base_polymer=polymer_generator, num_units=4, itp_path="polymer_n5.itp"
)
expander.expand_polymer(
    output_itp="polymer_n9.itp",
    output_pdb="polymer_n9.pdb",
    output_gro="polymer_n9.gro",
)
