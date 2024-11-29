from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from config.constants import ATOMISTIC_TRUNCATED_POLYMER_LENGTH
import os
from typing import Tuple


class AtomisticHomopolymerCreator:

    def __init__(self, num_units: float = ATOMISTIC_TRUNCATED_POLYMER_LENGTH):
        self.num_units = num_units

    @staticmethod
    def _create_polymerizable_monomer(monomer_smiles: str) -> Tuple[str, int, int]:
        """Create a polymerizable monomer by converting a double bond to single bonds."""
        monomer = Chem.MolFromSmiles(monomer_smiles)
        monomer = Chem.AddHs(monomer)  # Add explicit hydrogens

        # Identify the double bond to replace with single bonds
        double_bond = None
        for bond in monomer.GetBonds():
            if (
                bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                and not bond.GetIsAromatic()
                and not bond.IsInRing()
            ):
                double_bond = bond
                break

        if double_bond is None:
            raise ValueError(
                "No suitable non-aromatic double bond found for polymerization."
            )

        # Get the atom indices of the double bond
        end1_idx = double_bond.GetBeginAtomIdx()
        end2_idx = double_bond.GetEndAtomIdx()

        # Convert the double bond to a single bond to create open bonding sites
        polymerized_monomer = Chem.RWMol(monomer)
        polymerized_monomer.RemoveBond(end1_idx, end2_idx)
        polymerized_monomer.AddBond(end1_idx, end2_idx, Chem.rdchem.BondType.SINGLE)

        # Adjust hydrogens to ensure proper bonding sites
        polymerized_monomer.GetAtomWithIdx(end1_idx).SetNumExplicitHs(0)
        polymerized_monomer.GetAtomWithIdx(end2_idx).SetNumExplicitHs(0)
        Chem.SanitizeMol(polymerized_monomer)

        return polymerized_monomer, end1_idx, end2_idx

    def _create_linear_polymer(
        self, polymerized_monomer: str, end1_idx: float, end2_idx: int
    ):
        """Create a linear polymer by iteratively linking polymerizable monomers."""

        # Initialize the polymer with the first monomer
        polymer = Chem.RWMol(polymerized_monomer)

        # Track the end atom of the current chain to connect new monomers
        prev_end_idx = end1_idx  # Start with one end of the first monomer

        for _ in range(self.num_units - 1):
            # Create a new instance of the monomer
            new_monomer = Chem.RWMol(polymerized_monomer)

            # Map atoms of the new monomer to their new indices in the growing polymer
            atom_map = {}
            for atom in new_monomer.GetAtoms():
                new_idx = polymer.AddAtom(atom)
                atom_map[atom.GetIdx()] = new_idx

            # Add bonds within the new monomer based on its original connectivity
            for bond in new_monomer.GetBonds():
                begin_idx = atom_map[bond.GetBeginAtomIdx()]
                end_idx = atom_map[bond.GetEndAtomIdx()]
                bond_type = bond.GetBondType()
                polymer.AddBond(begin_idx, end_idx, bond_type)

            # Connect the previous monomer to the new monomer at the open ends
            polymer.AddBond(
                prev_end_idx, atom_map[end2_idx], Chem.rdchem.BondType.SINGLE
            )

            # Update the index of the new end for the next iteration
            prev_end_idx = atom_map[end1_idx]

        # Sanitize the polymer to finalize valency and structure
        Chem.SanitizeMol(polymer)

        return polymer

    @staticmethod
    def visualize_polymer(polymer):
        """Visualize the polymer molecule using RDKit without explicit hydrogens."""
        polymer_no_h = Chem.RemoveHs(polymer)  # Remove hydrogens for visualization
        return Draw.MolToImage(polymer_no_h, size=(800, 300))

    @staticmethod
    def _optimize_polymer(polymer):
        AllChem.EmbedMolecule(polymer)
        AllChem.UFFOptimizeMolecule(polymer)
        AllChem.AddHs(polymer)

    @staticmethod
    def save_polymer(polymer, polymer_name, output_dir):
        """Save the polymer molecule to a file in SDF format."""
        output_path = os.path.join(output_dir, f"{polymer_name}.pdb")
        with open(output_path, "w") as f:
            f.write(Chem.MolToPDBBlock(polymer))
        return output_path

    def create(self, monomer_smiles, polymer_name, output_path):
        """Create a linear polymer from a monomer SMILES string."""
        polymerized_monomer, end1_idx, end2_idx = self._create_polymerizable_monomer(
            monomer_smiles
        )
        polymer = self._create_linear_polymer(polymerized_monomer, end1_idx, end2_idx)
        self._optimize_polymer(polymer)
        self.save_polymer(polymer, polymer_name, output_path)
        return output_path
