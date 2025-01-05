from dataclasses import dataclass


@dataclass
class Solvent:
    """
    Represents a solvent for the simulation. Produced via .pdb
    """

    name: str
    molecular_weight: float
    density: float
    pdb_path: str  # Path to the solvent .pdb file
    pdb_molecule_name: str  # Name of the molecule in the .pdb file
    compressibility: str
    ### SEE IF YOU CAN USE PDB PARSER TO GET MOLECULE NAME
