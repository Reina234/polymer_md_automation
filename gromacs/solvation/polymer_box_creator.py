import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from typing import Optional, List
from config.constants import LengthUnits


class PolmerBoxResize:
    OUTPUT_GRO_NAME = "polymer_box.gro"
    UNITS = LengthUnits.NANOMETER

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def run(
        self,
        input_gro_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = [3.0, 3.0, 3.0],
        additional_notes: Optional[str] = None,
    ) -> str:
        solute_box_gro_path = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, self.OUTPUT_GRO_NAME
        )
        output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
        os.makedirs(output_dir, exist_ok=True)
        editconf_command = [
            "gmx",
            "editconf",
            "-f",
            input_gro_path,
            "-o",
            solute_box_gro_path,
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
                input_gro_path, solute_box_gro_path, box_size_nm, additional_notes
            )
        return solute_box_gro_path

    def metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        box_size: List[float],
        additional_notes=None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS - editconf",
            "defatils": f"created a polymer box of size {box_size}",
            "action(s)": f"used molecule at {input_file_path}, saved at {output_file_path}",
            "additional_notes": additional_notes,
        }

    def _update_metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        box_size: List[float],
        additional_notes: Optional[str] = None,
    ) -> None:
        metadata = self.metadata(
            input_file_path, output_file_path, box_size, additional_notes
        )
        self.metadata_tracker.add_step(step_name="GROMACS", details=metadata)
