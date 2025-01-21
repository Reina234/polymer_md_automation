from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple, List
from modules.shared.utils.file_utils import (
    check_file_does_not_exist,
    directory_exists_check_wrapper,
)
import os


class BasePolymerGenerator2(ABC):
    def __init__(self, cap_smiles: str = "[H]"):
        self.cap_smiles = cap_smiles

    def _create_monomer_residue(
        self, monomer_smiles: str
    ) -> Tuple[Chem.Mol, List[Tuple[int, int]]]:
        """
        Creates a monomer residue from the given SMILES string.

        :param monomer_smiles: The SMILES string representing the monomer.
        :return: A tuple containing the monomer residue and its open bonding sites.
        """
        monomer = Chem.MolFromSmiles(monomer_smiles)
        monomer = Chem.AddHs(monomer)

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
        open_sites = []

        bond = double_bonds[0]
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()

        rw_monomer.RemoveBond(atom1, atom2)
        rw_monomer.AddBond(atom1, atom2, Chem.rdchem.BondType.SINGLE)
        rw_monomer.GetAtomWithIdx(atom1).SetNumExplicitHs(0)
        rw_monomer.GetAtomWithIdx(atom2).SetNumExplicitHs(0)
        open_sites.append((atom1, atom2))

        Chem.SanitizeMol(rw_monomer)
        return rw_monomer, open_sites

    def _cap_termini(self, polymer: Chem.Mol, end1_idx: int, end2_idx: int) -> Chem.Mol:
        cap = Chem.MolFromSmiles(self.cap_smiles)
        rw_polymer = Chem.RWMol(polymer)

        cap_atom_idx = rw_polymer.AddAtom(cap.GetAtomWithIdx(0))
        rw_polymer.AddBond(end1_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)

        cap_atom_idx = rw_polymer.AddAtom(cap.GetAtomWithIdx(0))
        rw_polymer.AddBond(end2_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)

        Chem.SanitizeMol(rw_polymer)
        return rw_polymer

    def _finalise_molecule(self, mol: Chem.Mol, uff_optimise: bool = True) -> Chem.Mol:
        AllChem.EmbedMolecule(mol)
        if uff_optimise:
            AllChem.UFFOptimizeMolecule(mol)
        AllChem.AddHs(mol)
        return mol

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

    @abstractmethod
    def _generate_filename(self, **kwargs):
        pass

    @abstractmethod
    def _generate_polymer_rdkit(self, **kwargs) -> Chem.Mol:
        pass

    @abstractmethod
    def generate_polymer(self, **kwargs) -> str:
        pass


from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple, List, Set


class BasePolymerGenerator(ABC):
    def __init__(self):
        self.capping_atom = "[H]"  # Hydrogen cap (this will always be [H])
        self.monomer_residue_smiles: Set[str] = set()  # Store monomer SMILES
        self._possible_end_residues: Set[Chem.Mol] = (
            set()
        )  # Store end residue (hydrogen-added) molecules
        self.end_residue_smiles: Set[str] = None

    def _create_monomer_residue(
        self, monomer_smiles: str
    ) -> Tuple[Chem.Mol, List[int]]:
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
        open_sites = []

        bond = double_bonds[0]
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()

        rw_monomer.RemoveBond(atom1, atom2)
        rw_monomer.AddBond(atom1, atom2, Chem.rdchem.BondType.SINGLE)
        rw_monomer.GetAtomWithIdx(atom1).SetNumExplicitHs(0)
        rw_monomer.GetAtomWithIdx(atom2).SetNumExplicitHs(0)
        open_sites.append((atom1, atom2))

        Chem.SanitizeMol(rw_monomer)

        # Generate SMILES for the monomer
        monomer_smiles = Chem.MolToSmiles(rw_monomer, isomericSmiles=True)
        self.monomer_residue_smiles.add(monomer_smiles)

        # Generate end residues by adding hydrogen to each open site
        self._create_possible_end_residues(rw_monomer, open_sites)

        return rw_monomer, open_sites

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

    def _create_possible_end_residues(
        self, monomer: Chem.Mol, open_sites: List[Tuple[int, int]]
    ):
        """
        Adds hydrogen to each open site to create two possible end residues.
        :param monomer: The Chem.Mol object representing the monomer.
        :param open_sites: A list of open sites (atom indices).
        """
        # Create two possible end residues by adding hydrogen to each open site
        for site in open_sites:
            for atom_idx in site:
                # Add hydrogen to atom_idx
                modified_mol = Chem.RWMol(monomer)
                cap_atom_idx = modified_mol.AddAtom(Chem.Atom(1))  # Add hydrogen atom
                modified_mol.AddBond(
                    atom_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE
                )

                # Sanitize the molecule and add to end_caps
                Chem.SanitizeMol(modified_mol)
                self._possible_end_residues.add(modified_mol)

    def _cap_termini(self, polymer: Chem.Mol, end1_idx: int, end2_idx: int) -> Chem.Mol:
        cap = Chem.MolFromSmiles(self.capping_atom)
        rw_polymer = Chem.RWMol(polymer)

        cap_atom_idx = rw_polymer.AddAtom(cap.GetAtomWithIdx(0))
        rw_polymer.AddBond(end1_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)

        cap_atom_idx = rw_polymer.AddAtom(cap.GetAtomWithIdx(0))
        rw_polymer.AddBond(end2_idx, cap_atom_idx, Chem.rdchem.BondType.SINGLE)

        Chem.SanitizeMol(rw_polymer)
        return rw_polymer

    @staticmethod
    def residue_to_smiles(residue: Chem.Mol, open_sites: List[int]) -> str:
        """
        Converts a monomer residue into a SMILES string, adding placeholders for open sites.
        :param residue: The RDKit Mol object representing the monomer residue.
        :param open_sites: The open sites where polymerization can occur.
        :return: The SMILES string representing the residue.
        """
        rw_mol = Chem.RWMol(residue)

        for atom_idx in open_sites:
            placeholder_idx = rw_mol.AddAtom(
                Chem.Atom(0)
            )  # Adding placeholder atom (*)
            rw_mol.AddBond(
                atom_idx, placeholder_idx, Chem.rdchem.BondType.SINGLE
            )  # Bond placeholder to open site

        # Sanitize and convert the molecule to SMILES
        Chem.SanitizeMol(rw_mol)
        return Chem.MolToSmiles(rw_mol, isomericSmiles=True)

    def _match_end_residues(self, polymer: Chem.Mol) -> Set[str]:
        """
        Identifies the two end termini of a polymer by performing substructure matches with possible end residues.
        This function ensures we get exactly 2 matches, even if both matches correspond to the same end residue.
        :param polymer: The Chem.Mol object representing the polymer.
        :return: A set of matched end residue SMILES.
        """
        matched_end_residues = set()  # To store the unique matches
        match_count = 0  # To keep track of how many matches we've found

        # Iterate through all end residues and match with the polymer
        for end_residue in self._possible_end_residues:
            # Convert the end residue to SMILES for matching
            end_residue_smiles = Chem.MolToSmiles(end_residue, isomericSmiles=True)

            # Perform substructure matching
            matches = polymer.GetSubstructMatches(end_residue)

            if len(matches) > 0:
                # If we find matches, add them to the matched_termini set
                matched_end_residues.add(end_residue_smiles)
                match_count += len(matches)

        # Ensure exactly two termini are found (even if it's the same residue twice)
        if match_count != 2:
            raise ValueError(
                f"Expected exactly two end termini matches, but found {match_count}."
            )

        return matched_end_residues

    def _process_polymer(self, polymer: Chem.Mol):
        self.polymer_smiles = self.residue_to_smiles(polymer, [])
        matched_end_residues = self._match_end_residues(polymer)
        self.end_residue_smiles = matched_end_residues

    def retrieve_smiles(self):
        return self.monomer_residue_smiles, self.end_residue_smiles, self.polymer_smiles

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
