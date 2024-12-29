from abc import ABC, abstractmethod
from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Tuple, List
from A_modules.shared.utils.file_utils import (
    check_file_does_not_exist,
    directory_exists_check_wrapper,
)
import os


class BasePolymerGenerator(ABC):
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
