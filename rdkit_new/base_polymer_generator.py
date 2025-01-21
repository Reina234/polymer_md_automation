from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple, List, Set, Optional, Dict
import pandas as pd
from modules.shared.utils.file_utils import (
    check_file_does_not_exist,
    directory_exists_check_wrapper,
)
import os


class BasePolymerGenerator(ABC):
    def __init__(self):
        self.capping_atom = "[H]"  # Hydrogen cap (this will always be [H])
        self.map = []  # Mapping information as a list of dictionaries
        self.monomer_bead_map: Dict[str, str] = {}  # Map SMILES → Bead Type

    def _add_to_mapping(self, atom_indices: List[int], bead_type: str):
        """
        Adds atom indices to the mapping with the appropriate bead type.
        :param atom_indices: List of atom indices to map.
        :param bead_type: The bead type identifier.
        """
        current_count = sum(1 for entry in self.map if entry["bead_type"] == bead_type)
        for idx in atom_indices:
            unique_bead_name = f"{bead_type}{current_count + 1}"
            self.map.append(
                {
                    "atom_index": idx,
                    "bead_type": bead_type,
                    "unique_name": unique_bead_name,
                }
            )
            current_count += 1

    def _generate_bead_type(self, monomer_smiles: Optional[str] = None) -> str:
        """
        Generates a unique bead type for a given monomer SMILES.
        Ensures consistency when adding the same monomer multiple times.
        """
        if monomer_smiles is None or monomer_smiles not in self.monomer_bead_map:
            bead_count = len(self.monomer_bead_map) + 1
            self.monomer_bead_map[monomer_smiles] = f"B{bead_count}"
        return self.monomer_bead_map[monomer_smiles]

    def _create_monomer_residue(
        self, monomer_smiles: str
    ) -> Tuple[Chem.Mol, List[int], str]:
        """
        Creates a monomer residue from the given SMILES string, identifies open sites,
        and generates possible end residues with hydrogen.
        :param monomer_smiles: The SMILES string representing the monomer.
        :return: A tuple containing the monomer residue and its open bonding sites.
        """
        monomer = Chem.MolFromSmiles(monomer_smiles)
        monomer = Chem.AddHs(monomer)  # Add hydrogens to the molecule

        double_bonds = []
        for bond in monomer.GetBonds():
            if (
                bond.GetBondType() == Chem.rdchem.BondType.DOUBLE
                and not bond.GetIsAromatic()
                and not bond.IsInRing()
            ):
                double_bonds.append(bond)

        if len(double_bonds) > 1:
            raise ValueError("Monomer contains more than one double bond.")

        rw_monomer = Chem.RWMol(monomer)
        # open_sites = []

        bond = double_bonds[0]
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()

        rw_monomer.RemoveBond(atom1, atom2)
        rw_monomer.AddBond(atom1, atom2, Chem.rdchem.BondType.SINGLE)
        rw_monomer.GetAtomWithIdx(atom1).SetNumExplicitHs(0)
        rw_monomer.GetAtomWithIdx(atom2).SetNumExplicitHs(0)
        open_sites = (atom1, atom2)

        Chem.SanitizeMol(rw_monomer)

        bead_type = self._generate_bead_type(monomer_smiles)

        return rw_monomer, open_sites, bead_type

    def _add_monomer_to_polymer(
        self,
        polymer: Chem.RWMol,
        monomer: Chem.Mol,
        prev_end_idx: int,
        open_sites: Tuple[int, int],
        bead_type: str,
    ) -> Tuple[Chem.RWMol, int]:
        """
        Adds a monomer to the polymer at the specified open site.
        Ensures the polymer maintains exactly **one open site** after each addition.

        :param polymer: Growing polymer.
        :param monomer: Monomer residue (should have exactly **two open sites**).
        :param prev_end_idx: Index of the current open site in the growing polymer.
        :param open_sites: Tuple of **monomer's two open site indices**.
        :param bead_type: Bead type assigned to this monomer.
        :return: (Updated polymer, new open site index).
        """

        # Create a new copy of the monomer to add to the polymer
        new_monomer = Chem.RWMol(monomer)

        # Map atoms from new monomer into the growing polymer
        atom_map = {}
        for atom in new_monomer.GetAtoms():
            new_idx = polymer.AddAtom(atom)
            atom_map[atom.GetIdx()] = new_idx

        # Ensure **only one** open site remains after addition
        monomer_site_1 = atom_map[open_sites[0]]
        monomer_site_2 = atom_map[open_sites[1]]

        # Add intra-monomer bonds
        for bond in new_monomer.GetBonds():
            polymer.AddBond(
                atom_map[bond.GetBeginAtomIdx()],
                atom_map[bond.GetEndAtomIdx()],
                bond.GetBondType(),
            )

        # Attach the monomer to the polymer
        polymer.AddBond(prev_end_idx, monomer_site_1, Chem.rdchem.BondType.SINGLE)

        # Assign all atoms in this monomer to one bead
        new_indices = list(atom_map.values())
        self._add_to_mapping(new_indices, bead_type)

        # Return the updated polymer and **new open site index**
        return polymer, monomer_site_2

    def _create_cap_residues(
        self,
        monomer: Chem.Mol,
        open_sites: Tuple[int, int],
        use_open_site: int = 0,
        add_to_map: bool = True,
    ) -> Tuple[Chem.Mol, str, int]:
        """
        Caps the first monomer by adding a hydrogen (or specified capping atom) to the selected open site.

        :param monomer: The monomer molecule before capping.
        :param open_sites: Tuple of two open site atom indices.
        :param use_open_site: Index of the open site to cap (0 or 1).
        :param add_to_map: Whether to add this to the mapping.
        :return: (Hydrogen-capped monomer, bead type for mapping, new open site index).
        """
        if use_open_site not in [0, 1]:
            raise ValueError("Invalid open site index. Must be 0 or 1.")
        cap = Chem.MolFromSmiles(self.capping_atom)
        rw_monomer = Chem.RWMol(monomer)
        bead_type = self._generate_bead_type()

        # Choose which open site to cap
        cap_target_idx = open_sites[use_open_site]  # Ensure this is an integer
        unused_open_site = open_sites[1 - use_open_site]  # Get the other site

        cap_atom_idx = rw_monomer.AddAtom(cap.GetAtomWithIdx(0))
        rw_monomer.AddBond(cap_target_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)
        # Add the capping atom (default: hydrogen)

        Chem.SanitizeMol(rw_monomer)
        # Optionally add to the mapping
        if add_to_map:
            self._add_to_mapping([cap_target_idx, cap_atom_idx], f"{bead_type}_cap")

        return rw_monomer, bead_type, unused_open_site  # The remaining open site

    def _add_to_mapping(self, atom_indices: List[int], bead_type: str):
        """
        Adds atom indices to the mapping with the appropriate bead type.
        :param atom_indices: List of atom indices to map.
        :param bead_type: The bead type identifier.
        """
        current_count = sum(1 for entry in self.map if entry["bead_type"] == bead_type)
        for idx in atom_indices:
            unique_bead_name = f"{bead_type}{current_count + 1}"
            self.map.append(
                {
                    "atom_index": idx,
                    "bead_type": bead_type,
                    "unique_name": unique_bead_name,
                }
            )
            current_count += 1

    def _finalise_molecule(self, mol: Chem.Mol, uff_optimise: bool = True) -> Chem.Mol:
        """
        Finalizes the molecule by embedding it and adding hydrogens.
        :param mol: The molecule to finalize.
        :param uff_optimise: Flag to optimize the molecule using UFF.
        :return: The finalized molecule.
        """
        AllChem.EmbedMolecule(mol)
        if uff_optimise:
            AllChem.UFFOptimizeMolecule(mol)
        AllChem.AddHs(mol)
        return mol

    @abstractmethod
    def _generate_filename(self, **kwargs):
        pass

    @abstractmethod
    def _generate_polymer_rdkit(self, **kwargs) -> Chem.Mol:
        pass

    @abstractmethod
    def generate_polymer(self, **kwargs) -> str:
        pass

    @directory_exists_check_wrapper(dir_arg_index=2)
    def _save_as_pdb(
        self,
        polymer: Chem.Mol,
        output_dir: str,
        output_name: str,
        overwrite: bool = True,
    ) -> str:
        output_basename = output_name + ".pdb"
        output_path = os.path.join(output_dir, output_basename)
        if overwrite:
            suppress_error = True
            delete_file = True
        else:
            suppress_error = False
            delete_file = False
        check_file_does_not_exist(
            output_path, suppress_error=suppress_error, delete_file=delete_file
        )
        with open(output_path, "w") as f:
            f.write(Chem.MolToPDBBlock(polymer))
        return output_path

    @staticmethod
    def debug_print_mol(mol: Chem.Mol, label: str = "Molecule"):
        """
        Debugging function to print the SMILES of a molecule
        with and without explicit hydrogens, and check atom valencies.

        :param mol: The RDKit molecule to debug.
        :param label: A custom label for printing.
        """
        if not mol:
            print(f"[ERROR] {label}: Molecule is None!")
            return

        # Convert RWMol to Mol if needed
        if isinstance(mol, Chem.RWMol):
            mol = Chem.Mol(mol)

        mol_without_h = Chem.RemoveHs(mol)
        smiles_without_h = Chem.MolToSmiles(mol_without_h, isomericSmiles=True)

        # SMILES with explicit Hs
        mol_with_h = Chem.AddHs(mol)
        smiles_with_h = Chem.MolToSmiles(mol_with_h, isomericSmiles=True)

        print(
            f"[DEBUG] {label} (unchanged)   : {Chem.MolToSmiles(mol, isomericSmiles=True)}"
        )
        print(f"[DEBUG] {label} (No Hs)   : {smiles_without_h}")
        print(f"[DEBUG] {label} (With Hs) : {smiles_with_h}")

        # Print valency of each atom
        for atom in mol.GetAtoms():
            print(
                f"Atom {atom.GetIdx()} ({atom.GetSymbol()}): "
                f"Valency {atom.GetTotalValence()} | "
                f"Explicit Hs {atom.GetTotalNumHs()}"
            )
