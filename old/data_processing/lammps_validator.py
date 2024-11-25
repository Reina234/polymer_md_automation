import re
from typing import Dict


class LAMMPSValidator:
    """
    Validates the contents and format of LAMMPS topology and settings files.
    """

    @staticmethod
    def validate_pdb(pdb_file: str) -> None:
        """
        Validate the PDB file's existence and format.
        """
        with open(pdb_file, "r") as f:
            contents = f.read()
        if "ATOM" not in contents:
            raise ValueError(f"PDB file {pdb_file} is invalid or corrupted.")
        print(f"PDB validation passed: {pdb_file}")

    @staticmethod
    def validate_topology(topology_files: Dict[str, str]) -> None:
        """
        Validate the LAMMPS data and settings files.
        """
        required_sections = ["Atoms", "Bonds"]
        for key, file_path in topology_files.items():
            with open(file_path, "r") as f:
                contents = f.read()
            for section in required_sections:
                if not re.search(f"^\\s*{section}\\s*$", contents, re.MULTILINE):
                    raise ValueError(f"Section '{section}' is missing in {key} file: {file_path}")
        print("Topology validation passed.")
