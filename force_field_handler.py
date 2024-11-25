import os
import subprocess
class ForceFieldHandler:
    """
    Handles parameterization of molecules.
    """

    def __init__(self, method: str = "GAFF"):
        self.method = method

    def parameterize_with_acpype(self, pdb_file: str, output_dir: str) -> str:
        """
        Use ACPYPE to parameterize the molecule with GAFF.
        """
        acpype_output = os.path.join(output_dir, f"{os.path.basename(pdb_file).split('.')[0]}_acpype")
        subprocess.run(
            [
                "acpype",
                "-i", pdb_file,
                "-o", acpype_output,
                "--gmx",
            ],
            check=True
        )
        return acpype_output
