import os
import subprocess
from processing.metadata_tracker import MetadataTracker
from config.gromacs_paths import GROMACS_ION_SCRIPT


class IonAdder:
    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def add_ions(self, input_gro: str, topology_file: str, output_dir: str, mdp_file: str = GROMACS_ION_SCRIPT) -> dict:
        """
        Add ions to neutralize the simulation box.

        Args:
            input_gro (str): Path to the solvated .gro file.
            topology_file (str): Path to the topology file.
            mdp_file (str): Path to the ions.mdp file.
            output_dir (str): Directory to save output files.

        Returns:
            dict: Path to the neutralized .gro file and updated topology file.
        """
        os.makedirs(output_dir, exist_ok=True)
        tpr_file = os.path.join(output_dir, "ions.tpr")
        output_gro = os.path.join(output_dir, "polymer_neutralized.gro")

        # Generate .tpr file
        grompp_command = [
            "gmx", "grompp",
            "-f", mdp_file,
            "-c", input_gro,
            "-p", topology_file,
            "-o", tpr_file
        ]
        subprocess.run(grompp_command, check=True)

        # Add ions
        genion_command = [
            "gmx", "genion",
            "-s", tpr_file,
            "-o", output_gro,
            "-p", topology_file,
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ]
        subprocess.run(genion_command, check=True)

        self.metadata_tracker.add_step("Ion Addition", {
            "input_gro": input_gro,
            "topology_file": topology_file,
            "output_gro": output_gro
        })

        return {"neutralized_gro_file": output_gro, "topology_file": topology_file}
