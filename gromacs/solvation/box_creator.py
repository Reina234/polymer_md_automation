import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_DIR
from typing import Optional, Tuple


class BoxCreator:
    OUTPUT_GRO_NAME = "polymer_box.gro"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def create_box(
        self,
        input_gro_path: str,
        output_dir: str = GROMACS_OUTPUT_DIR,
        box_size_nm: tuple = (2.0, 2.0, 2.0),
        additional_notes: Optional[str] = None,
    ) -> str:
        output_file_path = os.path.join(output_dir, self.OUTPUT_GRO_NAME)
        os.makedirs(output_dir, exist_ok=True)
        editconf_command = [
            "gmx",
            "editconf",
            "-f",
            input_gro_path,
            "-o",
            output_file_path,
            "-box",
            str(box_size_nm[0]),
            str(box_size_nm[1]),
            str(box_size_nm[2]),
            "-c",  # Distance to the box edge
            "-bt",
            "cubic",  # Box type
        ]
        subprocess.run(editconf_command, check=True)
        if self.metadata_tracker:
            self._update_metadata(
                input_gro_path, output_file_path, box_size_nm, additional_notes
            )
        return output_file_path

    def metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        box_size: Tuple[str],
        additional_notes=None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS",
            "defatils": f"created a solvent box of size {box_size}",
            "action(s)": f"used molecule at {input_file_path}, saved at {output_file_path}",
            "additional_notes": additional_notes,
        }

    def _update_metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        box_size: Tuple[str],
        additional_notes: Optional[str] = None,
    ) -> None:
        metadata = self.metadata(
            input_file_path, output_file_path, box_size, additional_notes
        )
        self.metadata_tracker.add_step(step_name="GROMACS", details=metadata)
