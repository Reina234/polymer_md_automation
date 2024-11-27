import os
import subprocess
from processing.metadata_tracker import MetadataTracker
from data_models.solvent import Solvent


class SolvateBox:
    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def solvate(self, box_file: str, solvent: Solvent, polymer_top_file: str, output_dir: str) -> dict:
        """
        Add solvent to the box using GROMACS solvate.

        Args:
            box_file (str): Path to the cubic box file.
            solvent (Solvent): Solvent object.
            polymer_top_file (str): Path to the polymer topology file.
            output_dir (str): Directory to save output files.

        Returns:
            dict: Paths to the solvated box and updated topology file.
        """
        solvated_gro_file = os.path.join(output_dir, "polymer_solvated.gro")

        solvate_command = [
            "gmx", "solvate",
            "-cp", box_file,
            "-cs", solvent.pdb_path,
            "-o", solvated_gro_file,
            "-p", polymer_top_file
        ]
        subprocess.run(solvate_command, check=True)

        self.metadata_tracker.add_step("Solvation", {
            "box_file": box_file,
            "solvent_pdb": solvent.pdb_path,
            "solvated_gro_file": solvated_gro_file,
            "topology_file": polymer_top_file
        })

        return {"solvated_gro_file": solvated_gro_file, "topology_file": polymer_top_file}
