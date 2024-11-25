from dataclasses import dataclass

@dataclass(frozen=True)
class Solvent:
    """
    Represents a solvent
    """
    name: str
    molecular_weight: float
    density: float

