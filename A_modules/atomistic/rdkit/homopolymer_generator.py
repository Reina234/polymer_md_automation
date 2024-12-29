from rdkit import Chem
from typing import Optional
from A_modules.atomistic.rdkit.base_polymer_generator import BasePolymerGenerator


class HomopolymerGenerator(BasePolymerGenerator):
    def __init__(self, cap_smiles: str = "[H]"):
        super().__init__(cap_smiles)

    def _generate_polymer_rdkit(self, monomer_smiles: str, num_units: int) -> Chem.Mol:
        monomer_residue, open_sites = self._create_monomer_residue(monomer_smiles)

        polymer = Chem.RWMol(monomer_residue)

        prev_end_idx = open_sites[0][0]
        last_end_idx = open_sites[0][1]

        for _ in range(num_units - 1):
            new_monomer = Chem.RWMol(monomer_residue)

            atom_map = {}
            for atom in new_monomer.GetAtoms():
                new_idx = polymer.AddAtom(atom)
                atom_map[atom.GetIdx()] = new_idx

            for bond in new_monomer.GetBonds():
                begin_idx = atom_map[bond.GetBeginAtomIdx()]
                end_idx = atom_map[bond.GetEndAtomIdx()]
                bond_type = bond.GetBondType()
                polymer.AddBond(begin_idx, end_idx, bond_type)

            polymer.AddBond(
                prev_end_idx, atom_map[open_sites[0][1]], Chem.rdchem.BondType.SINGLE
            )
            prev_end_idx = atom_map[open_sites[0][0]]

        polymer = self._cap_termini(
            polymer=polymer, end1_idx=prev_end_idx, end2_idx=last_end_idx
        )
        return polymer

    def _generate_filename(self, monomer_smiles: str, num_units: int) -> str:
        return f"{monomer_smiles.lower()}_{num_units}"

    def generate_polymer(
        self,
        monomer_smiles: str,
        num_units: int,
        output_dir: str,
        output_name: Optional[str] = None,
        uff_optimise: bool = True,
        overwrite: bool = True,
    ) -> str:
        polymer = self._generate_polymer_rdkit(monomer_smiles, num_units)
        polymer = self._finalise_molecule(polymer, uff_optimise=uff_optimise)
        if output_name is None:
            output_name = self._generate_filename(monomer_smiles, num_units)
        else:
            output_name = f"{output_name}.pdb"
        output_path = self._save_as_pdb(
            polymer,
            output_dir,
            output_name=output_name,
            overwrite=overwrite,
        )
        return output_path
