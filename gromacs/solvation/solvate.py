import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from typing import Optional, List
from preprocessing.pdb_utils import calculate_num_particles
from config.constants import LengthUnits
from pathlib import Path
from preprocessing.utils import move_and_rename_topol_file

# NOTE: add in file type checks


class Solvate:
    OUTPUT_GRO_NAME = "solvated_polymer.gro"
    TOPOL_NAME = "topol.top"
    UNITS = LengthUnits.NANOMETER

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def run(
        self,
        solvent_box_gro_path: str,
        solute_box_gro_path: str,
        input_top_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        additional_notes: Optional[str] = None,
    ):

        output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
        os.makedirs(output_dir, exist_ok=True)
        topol_file = move_and_rename_topol_file(
            input_top_path, output_dir, self.TOPOL_NAME
        )
        solvated_solute_gro_path = os.path.join(output_dir, self.OUTPUT_GRO_NAME)
        solvate_command = [
            "gmx",
            "solvate",
            "-cp",
            solute_box_gro_path,
            "-cs",
            solvent_box_gro_path,
            "-o",
            solvated_solute_gro_path,
            "-p",
            topol_file,
        ]
        subprocess.run(solvate_command, check=True)
        if self.metadata_tracker:
            self._update_metadata(
                solvent_box_gro_path,
                solute_box_gro_path,
                solvated_solute_gro_path,
                additional_notes,
            )
        return solvated_solute_gro_path

    def metadata(
        self,
        solvent_box_gro_path: str,
        solute_box_gro_path: str,
        output_file_path: str,
        additional_notes=None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS solvate",
            "details": f"added solute at {solute_box_gro_path} to {solvent_box_gro_path}",
            "action(s)": f"aved at {output_file_path}",
            "additional_notes": additional_notes,
        }

    def _update_metadata(
        self,
        solvent_box_gro_path: str,
        solute_box_gro_path: str,
        output_file_path: str,
        additional_notes=None,
    ) -> None:
        metadata = self.metadata(
            solvent_box_gro_path,
            solute_box_gro_path,
            output_file_path,
            additional_notes,
        )
        self.metadata_tracker.add_step(step_name="GROMACS", details=metadata)
