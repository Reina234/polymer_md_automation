from dataclasses import dataclass

@dataclass
class Solvent:
    """
    Represents a solvent for the simulation.
    """
    name: str
    molecular_weight: float
    density: float
    itp_path: str  # Path to the solvent .itp file
    pdb_path: str  # Path to the solvent .pdb file
    force_field: str  # Path to the force field directory
