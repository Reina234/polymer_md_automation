# providers/topology_moltemplate_provider.py
import os
import subprocess
from .base_topology_provider import TopologyProvider


class TopologyFromMoltemplate(TopologyProvider):
    def __init__(self, lt_file: str, output_dir: str = "moltemplate_files"):
        self.lt_file = lt_file
        self.output_dir = output_dir

    def get_topology(self) -> Dict:
        """
        Generate topology and parameters from a Moltemplate file.
        """
        os.makedirs(self.output_dir, exist_ok=True)
        result = subprocess.run(
            ["moltemplate.sh", self.lt_file],
            cwd=self.output_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Moltemplate failed: {result.stderr}")

        data_file = os.path.join(self.output_dir, f"{os.path.splitext(os.path.basename(self.lt_file))[0]}.data")
        settings_file = os.path.join(self.output_dir, f"{os.path.splitext(os.path.basename(self.lt_file))[0]}.settings")
        return {"data_file": data_file, "settings_file": settings_file}
