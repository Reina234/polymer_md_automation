from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem import AllChem
from typing import Tuple, List, Set, Optional, Dict
import pandas as pd
from modules.shared.utils.file_utils import (
    check_file_does_not_exist,
    directory_exists_check_wrapper,
)
import os
from itertools import cycle


class BasePolymerGenerator(ABC):

    def __init__(self, monomer_smiles: List[str]):
        self.capping_atom = "[H]"  # Hydrogen cap (this will always be [H])
        self.pycg_map = []  # Mapping information as a list of dictionaries
        self.votca_map = []  # Mapping information for VOTCA
        self.votca_bonds = []  # List of VOTCA bond tuples (prev_bead, new_bead)
        self.votca_angles = []  # List of angle triplets for VOTCA
        self.monomer_bead_map: Dict[str, str] = {}  # Mapping SMILES to bead types
        self.bead_type_counts: Dict[str, int] = {}  # Tracking counts per bead type
        self.mol: Optional[Chem.Mol] = None
        self.monomer_smiles_list: List[str] = monomer_smiles
        self.sequence: List[str] = []
        self.pdb_path: Optional[str] = None

    def _add_to_sequence(self, mol: Chem.Mol) -> None:
        mol_smiles = Chem.MolToSmiles(mol)
        self.sequence.append(mol_smiles)

    def _generate_bead_type(self, monomer_smiles: Optional[str] = None) -> str:
        """
        Generates a bead type for a given monomer SMILES.
        Ensures consistency when adding the same monomer multiple times.
        """
        if monomer_smiles not in self.monomer_bead_map:
            bead_count = len(self.monomer_bead_map) + 1
            self.monomer_bead_map[monomer_smiles] = f"B{bead_count}"

        return self.monomer_bead_map[monomer_smiles]

    def _generate_unique_bead_name(self, bead_type: str) -> str:
        """
        Generates a unique bead name for a given bead type.
        Beads are named as 'B1N1', 'B1N2', 'B2N1', etc.

        :param bead_type: The base bead type identifier.
        :return: A unique bead name.
        """
        if bead_type not in self.bead_type_counts:
            self.bead_type_counts[bead_type] = 1
        else:
            self.bead_type_counts[bead_type] += 1

        return f"{bead_type}N{self.bead_type_counts[bead_type]}"

    def _add_to_mapping(
        self, atom_indices: List[int], bead_type: str, unique_name: str, mol: Chem.Mol
    ):
        """
        Adds atom indices to both self.cg_map (for PyCGTOOL) and self.votca_map (for VOTCA).
        Ensures that PyCGTOOL mapping is stored per atom, while VOTCA mapping is stored per bead.

        :param atom_indices: List of atom indices belonging to a monomer.
        :param bead_type: The bead type identifier.
        :param unique_name: The unique bead name for this monomer residue.
        """
        self._add_to_pycg_map(atom_indices, bead_type, unique_name)
        self._add_to_votca_map(atom_indices, bead_type, unique_name, mol)

    def _add_bond_to_votca(self, new_bead: str):
        """
        Adds a bond to the VOTCA bond list.
        Ensures that the new bead bonds with the last added bead.
        """
        if len(self.votca_map) < 1:
            return  # No previous bead to bond with

        prev_bead = self.votca_map[-1]["unique_name"]  # Last added bead before this one
        if (prev_bead, new_bead) not in self.votca_bonds:
            self.votca_bonds.append((prev_bead, new_bead))

    def _add_angle_to_votca(self, new_bead: str):
        """
        Adds an angle to the VOTCA angle list.
        Uses the last three beads added to define an angle.
        """
        if len(self.votca_bonds) > 1:
            bead1, bead2 = self.votca_bonds[-1]  # First bond
            bead3 = new_bead
            self.votca_angles.append((bead1, bead2, bead3))

    def _add_to_pycg_map(
        self, atom_indices: List[int], bead_type: str, unique_name: str
    ):
        """
        Stores mapping data for PyCGTOOL in self.cg_map.
        Each atom is recorded individually.

        :param atom_indices: List of atom indices in a monomer.
        :param bead_type: The bead type identifier.
        :param unique_name: Unique name assigned to this bead.
        """
        for idx in atom_indices:
            self.pycg_map.append(
                {
                    "atom_index": idx,  # Each atom has a separate row
                    "bead_type": bead_type,
                    "unique_name": unique_name,
                }
            )

    def _add_to_votca_map(
        self, atom_indices: List[int], bead_type: str, unique_name: str, mol: Chem.Mol
    ):
        """
        Stores mapping data for VOTCA, including atom indices and predicted atom names.

        :param atom_indices: List of atom indices in a monomer.
        :param bead_type: The bead type identifier.
        :param unique_name: Unique bead name assigned to this monomer residue.
        """
        mol_smiles = Chem.MolToSmiles(mol)
        if not atom_indices or not self.mol:
            return

        periodic_table = Chem.GetPeriodicTable()

        # Compute atom names & masses
        atom_names = []
        atomic_masses = []
        for idx in atom_indices:
            atom = self.mol.GetAtomWithIdx(idx)
            atom_names.append(atom.GetSymbol() + str(idx))  # Predict atom name
            atomic_masses.append(periodic_table.GetAtomicWeight(atom.GetSymbol()))

        # Compute mass fractions (normalized)

        mass_weights = [round(mass) for mass in atomic_masses]

        # Store mapping
        self.votca_map.append(
            {
                "unique_name": unique_name,
                "bead_type": bead_type,
                "atom_indices": atom_indices,  # Atom indices
                "atom_names": atom_names,  # Predicted atom names
                "x-weight": mass_weights,  # Mass-weighted position
                "f-weight": mass_weights,  # Mass-weighted force distribution
                "smiles": mol_smiles,
            }
        )

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
        add_to_sequence: bool = True,
    ) -> Tuple[Chem.RWMol, int]:
        """
        Adds a monomer to the polymer at the specified open site.
        Ensures that the polymer maintains exactly **one open site** after each addition.

        :param polymer: Growing polymer.
        :param monomer: Monomer residue (should have exactly **two open sites**).
        :param prev_end_idx: Index of the current open site in the growing polymer.
        :param open_sites: Tuple of **monomer's two open site indices**.
        :param bead_type: Bead type assigned to this monomer.
        :return: (Updated polymer, new open site index).
        """

        # Create a new copy of the monomer to add to the polymer
        new_monomer = Chem.RWMol(monomer)
        if add_to_sequence:
            self._add_to_sequence(monomer)
        # Generate a unique name for this monomer residue
        unique_name = self._generate_unique_bead_name(bead_type)

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

        # Assign all atoms in this monomer to the same unique bead
        new_indices = list(atom_map.values())

        self.mol = polymer
        new_bead = unique_name
        self._add_bond_to_votca(new_bead)
        self._add_angle_to_votca(new_bead)  # Add angle to previous beads

        self._add_to_mapping(new_indices, bead_type, unique_name, monomer)

        # Return the updated polymer and **new open site index**
        return polymer, monomer_site_2

    def _create_cap_residues(
        self,
        monomer: Chem.Mol,
        open_sites: Tuple[int, int],
        use_open_site: int = 0,
        add_to_map: bool = True,
        add_to_sequence: bool = False,
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

        # Convert monomer to RWMol for modification
        rw_monomer = Chem.RWMol(monomer)

        # Generate a bead type and unique name for this capped residue

        # Determine which open site to cap and which remains open
        cap_target_idx = open_sites[use_open_site]  # The site to be capped
        unused_open_site = open_sites[1 - use_open_site]  # The site that remains open

        # Create and add the capping atom (default: hydrogen)
        cap_atom = Chem.Atom(1)  # Hydrogen
        cap_atom_idx = rw_monomer.AddAtom(cap_atom)

        # Bond the capping atom to the selected open site
        rw_monomer.AddBond(cap_target_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)

        # **Ensure all atoms in this monomer are mapped to a unique bead**
        all_atom_indices = [atom.GetIdx() for atom in rw_monomer.GetAtoms()]
        self.mol = rw_monomer
        cap_smiles = Chem.MolToSmiles(rw_monomer)
        bead_type = self._generate_bead_type(cap_smiles)
        unique_name = self._generate_unique_bead_name(bead_type)

        # Perform controlled sanitization, avoiding unintended changes
        Chem.SanitizeMol(
            rw_monomer,
            sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
            ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
        )
        if add_to_map:
            self._add_to_mapping(all_atom_indices, bead_type, unique_name, rw_monomer)

        if add_to_sequence:
            self._add_to_sequence(rw_monomer)
        return rw_monomer, bead_type, unused_open_site  # The remaining open site

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
    def _generate_polymer_rdkit(self, num_units: int) -> Chem.Mol:
        pass

    def generate_polymer(
        self,
        num_units: int,
        output_dir: str,
        output_name: Optional[str] = None,
        uff_optimise: bool = True,
        overwrite: bool = True,
        save: bool = True,
    ) -> str:
        polymer = self._generate_polymer_rdkit(num_units)
        polymer = self._finalise_molecule(polymer, uff_optimise=uff_optimise)
        if save:
            if output_name is None:
                output_name = self._generate_filename(num_units)
            else:
                output_name = f"{output_name}.pdb"
            output_path = self._save_as_pdb(
                polymer,
                output_dir,
                output_name=output_name,
                overwrite=overwrite,
            )
            return output_path
        else:
            return polymer

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

        self.pdb_path = output_path
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
