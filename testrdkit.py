from rdkit import Chem
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from modules.atomistic.rdkit.homopolymer_generator import HomopolymerGenerator
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple


def itp_to_dataframe(itp_file: str):
    parser = GromacsParser()
    sections = parser.parse("rdkit_test/POLY_GMX.itp")
    atom_section = sections["data_atoms"]
    data_handler = DataHandler()
    # data_handler.section = atom_section
    data_handler.process(atom_section)
    return data_handler.content


def match_atoms_to_beads(
    itp_df: pd.DataFrame, smiles_list: List[str]
) -> Dict[int, str]:
    """
    Match atoms from the .itp DataFrame to beads based on SMILES strings,
    considering explicit hydrogen counts.
    Args:
        itp_df (pd.DataFrame): Parsed .itp file as a DataFrame.
        smiles_list (List[str]): List of SMILES strings for the beads.
    Returns:
        Dict[int, str]: Mapping of atom indices to bead names.
    """
    atom_to_bead = {}
    monomer_mols = {
        f"BEAD_{i+1}": Chem.MolFromSmiles(smiles)
        for i, smiles in enumerate(smiles_list)
    }

    for bead_name, bead_mol in monomer_mols.items():
        if bead_mol is None:
            raise ValueError(f"Invalid SMILES for bead: {bead_name}")

        # Match atoms in the .itp to this bead, prioritizing end monomers
        for idx, row in itp_df.iterrows():
            atom_name = row["atom"]
            atom_type = row["type"]
            atom_h_count = row.get(
                "hydrogens", 0
            )  # Ensure explicit hydrogens are counted

            # Verify bead matches by substructure and hydrogen count
            if (
                bead_mol.HasSubstructMatch(Chem.MolFromSmiles(atom_name))
                and bead_mol.GetNumAtoms() == len(atom_name) + atom_h_count
            ):
                if idx + 1 not in atom_to_bead:
                    atom_to_bead[idx + 1] = bead_name

    return atom_to_bead


def generate_mapping_file(
    atom_to_bead: Dict[int, str], itp_df: pd.DataFrame, output_file: str = "mapping.txt"
):
    """
    Generate the PyCGTOOL mapping file.
    Args:
        atom_to_bead (Dict[int, str]): Mapping of atom indices to bead names.
        itp_df (pd.DataFrame): Parsed .itp DataFrame for reference.
        output_file (str): Path to the output mapping file.
    """
    with open(output_file, "w") as f:
        f.write("[MAPPING]\n")
        assigned_atoms = set()
        for atom_idx, bead_name in atom_to_bead.items():
            atom_name = itp_df.loc[atom_idx - 1, "atom"]  # Adjust index for DataFrame
            f.write(f"{bead_name} P3 {atom_name}\n")
            assigned_atoms.add(atom_idx)

        # Ensure all atoms are assigned
        unassigned_atoms = set(itp_df.index + 1) - assigned_atoms
        if unassigned_atoms:
            raise ValueError(f"Unassigned atoms detected: {unassigned_atoms}")

    print(f"Mapping file '{output_file}' generated successfully.")


def preprocess_smiles(smiles_list: List[str]) -> List[str]:
    """
    Preprocess and order SMILES strings for robust matching.
    Args:
        smiles_list (List[str]): List of SMILES strings.
    Returns:
        List[str]: Ordered list of SMILES strings.
    """
    return sorted(
        smiles_list,
        key=lambda x: (
            Chem.MolFromSmiles(x).GetNumAtoms(),  # Number of atoms
            len(x),  # String length as a proxy for complexity
        ),
        reverse=True,
    )


def verify_assignments(atom_to_bead: Dict[int, str], itp_df: pd.DataFrame) -> None:
    """
    Verify that all atoms are assigned to a bead.
    Args:
        atom_to_bead (Dict[int, str]): Mapping of atom indices to bead names.
        itp_df (pd.DataFrame): Parsed .itp DataFrame for reference.
    """
    all_atoms = set(itp_df.index + 1)
    assigned_atoms = set(atom_to_bead.keys())
    unassigned_atoms = all_atoms - assigned_atoms

    if unassigned_atoms:
        raise ValueError(f"Unassigned atoms detected: {unassigned_atoms}")


def hierarchical_match(
    polymer_mol: Chem.Mol, smiles_list: List[str], itp_df: pd.DataFrame
) -> Dict[int, str]:
    """
    Match atoms in the polymer to SMILES strings hierarchically.
    Args:
        polymer_mol (Chem.Mol): RDKit molecule of the entire polymer.
        smiles_list (List[str]): List of SMILES strings.
        itp_df (pd.DataFrame): Parsed .itp DataFrame for reference.
    Returns:
        Dict[int, str]: Mapping of atom indices (1-based) to molecule names (e.g., UNL).
    """
    atom_to_bead = {}
    assigned_atoms = set()

    # Preprocess and order SMILES
    ordered_smiles = preprocess_smiles(smiles_list)

    # Match in hierarchical order
    for smiles in ordered_smiles:
        monomer_mol = Chem.MolFromSmiles(smiles)
        if monomer_mol is None:
            raise ValueError(f"Invalid SMILES: {smiles}")

        matches = polymer_mol.GetSubstructMatches(monomer_mol)
        for match in matches:
            match_indices = [idx + 1 for idx in match]  # Convert to 1-based indices

            # Check for overlap
            overlap = any(idx in assigned_atoms for idx in match_indices)
            if not overlap:  # Only assign non-overlapping matches
                for atom_idx in match_indices:
                    atom_to_bead[atom_idx] = smiles
                    assigned_atoms.add(atom_idx)

    # Verify all atoms are assigned
    verify_assignments(atom_to_bead, itp_df)

    return atom_to_bead


generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 5, "rdkit_test2", overwrite=False, save=False)
smiles_list = generator.retrieve_unit_smiles_list()

itp_df = itp_to_dataframe("rdkit_test/POLY_GMX.itp")
# atom_to_bead = match_atoms_to_beads(itp_df, smiles_list)
atom_to_bead = hierarchical_match(generator.polymer_mol, smiles_list, itp_df)
print(atom_to_bead)
generate_mapping_file(atom_to_bead, itp_df, "map_test.txt")
