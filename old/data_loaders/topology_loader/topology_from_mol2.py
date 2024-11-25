import parmed as pmd
from .base_topology_provider import TopologyProvider
import os
from typing import Dict


class TopologyFromMOL2(TopologyProvider):
    def __init__(self, mol2_file: str):
        self.mol2_file = mol2_file

    def get_topology(self, output_dir: str) -> Dict[str, str]:
        """
        Generate LAMMPS topology files from a MOL2 file.
        """
        os.makedirs(output_dir, exist_ok=True)
        data_file = os.path.join(output_dir, "system.data")
        settings_file = os.path.join(output_dir, "system.settings")

        structure = pmd.load_file(self.mol2_file)
        structure.save(data_file, format="LAMMPSDATA")
        structure.save(settings_file, format="LAMMPS")
        print(f"LAMMPS topology files generated: {data_file}, {settings_file}")

        return {"data_file": data_file, "settings_file": settings_file}
