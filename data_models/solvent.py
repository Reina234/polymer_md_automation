from dataclasses import dataclass


@dataclass
class Solvent:
    """
    Represents a solvent for the simulation.
    """

    name: str
    molecular_weight: float
    density: float
    pdb_path: str  # Path to the solvent .pdb file
    pdb_molecule_name: str  # Name of the molecule in the .pdb file
