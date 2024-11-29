import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from typing import Optional, List
from config.constants import LengthUnits
from gromacs.base_gromacs_command import BaseGromacsCommand


class PolmerBoxResize(BaseGromacsCommand):
    OUTPUT_NAME = "polymer_box.gro"

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
        command, output_gro_path = self._create_grommp_command(
            input_gro_path=input_gro_path,
            run_name=run_name,
            output_base_dir=output_base_dir,
            box_size_nm=box_size_nm,
            additional_notes=additional_notes,
        )
        self._execute(command)
        if self.metadata_tracker:
            self._update_metadata(
                input_gro_path=input_gro_path,
                run_name=run_name,
                output_base_dir=output_base_dir,
                box_size_nm=box_size_nm,
                additional_notes=additional_notes,
            )

        return output_gro_path

    def _create_grommp_command(
        self,
        input_gro_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = [4.0, 4.0, 4.0],
        additional_notes: Optional[str] = None,
    ) -> str:
        solute_box_gro_path = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, self.OUTPUT_NAME
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
        return editconf_command, solute_box_gro_path

    def metadata(
        self,
        input_gro_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = [3.0, 3.0, 3.0],
        additional_notes: Optional[str] = None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS - editconf",
            "defatils": f"created a polymer box of size {box_size_nm} with units {self.UNITS.value}",
            "action(s)": f"used molecule at {input_gro_path}, saved at {output_base_dir}/{run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.OUTPUT_NAME}",
            "additional_notes": additional_notes,
        }
