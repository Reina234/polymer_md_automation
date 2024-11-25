# core/solvent.py
from dataclasses import dataclass
from core.pdb_provider import PDBProvider
from core.topology_provider import TopologyProvider


@dataclass
class Solvent:
    name: str
    molecular_weight: float
    density: float
    pdb_provider: PDBProvider
    topology_provider: TopologyProvider

    def generate_pdb(self, output_file: str) -> None:
        """
        Generate or load the PDB file using the assigned provider.
        """
        self.pdb_provider.get_pdb(output_file)

    def get_topology(self) -> Dict:
        """
        Generate or load the topology using the assigned provider.
        """
        return self.topology_provider.get_topology()
