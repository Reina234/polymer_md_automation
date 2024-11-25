import os
import subprocess
from typing import Dict


class MoltemplateAdapter(LAMMPSDataAdapter):
    """
    Adapter for converting Moltemplate files to LAMMPS-compatible data.
    """

    def __init__(self, lt_file: str):
        self.lt_file = lt_file

    def load(self, output_dir: str) -> Dict[str, str]:
        """
        Generate LAMMPS data and settings files from a Moltemplate file.
        """
        os.makedirs(output_dir, exist_ok=True)
        result = subprocess.run(
            ["moltemplate.sh", self.lt_file],
            cwd=output_dir,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(f"Moltemplate conversion failed:\n{result.stderr}")

        data_file = os.path.join(output_dir, "system.data")
        settings_file = os.path.join(output_dir, "system.settings")
        print(f"Loaded Moltemplate file. Data: {data_file}, Settings: {settings_file}")

        return {"data_file": data_file, "settings_file": settings_file}
