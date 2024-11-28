import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker


class BoxCreator:
    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def create_box(self, input_gro: str, output_dir: str, box_size: float) -> str:
        """
        Create a cubic box using GROMACS editconf.

        Args:
            input_gro (str): Path to the input .gro file.
            output_dir (str): Directory to save the output.
            box_size (float): Desired box size in nm.

        Returns:
            str: Path to the created box file.
        """
        box_file = os.path.join(output_dir, "polymer_box.gro")
        os.makedirs(output_dir, exist_ok=True)

        editconf_command = [
            "gmx",
            "editconf",
            "-f",
            input_gro,
            "-o",
            box_file,
            "-c",  # Center the polymer
            "-d",
            str(box_size),  # Distance to the box edge
            "-bt",
            "cubic",  # Box type
        ]
        subprocess.run(editconf_command, check=True)

        self.metadata_tracker.add_step(
            "Box Creation", {"box_file": box_file, "box_size": box_size}
        )
        return box_file
