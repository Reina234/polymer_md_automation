import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from typing import Optional, List
from preprocessing.calculation_utils import calculate_num_particles
from config.constants import LengthUnits
from pathlib import Path
from gromacs.gromacs_utils import move_and_rename_topol_file
from gromacs.base_gromacs_command import BaseGromacsCommand

# NOTE: add in file type checks


class Solvate(BaseGromacsCommand):
    OUTPUT_GRO_NAME = "solvated_polymer.gro"
    TOPOL_NAME = "topol.top"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def run(
        self,
        solvent_box_gro_path: str,
        solute_box_gro_path: str,
        input_top_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        additional_notes: Optional[str] = None,
    ):
        """
        Solvates a solute in a solvent box using GROMACS solvate.

        """

        return self._execute(
            solvent_box_gro_path=solvent_box_gro_path,
            solute_box_gro_path=solute_box_gro_path,
            input_top_path=input_top_path,
            run_name=run_name,
            output_base_dir=output_base_dir,
            additional_notes=additional_notes,
        )

    def _create_subprocess_command(
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
            input_top_path,
        ]
        return solvate_command, solvated_solute_gro_path

    def metadata(
        self,
        solvent_box_gro_path: str,
        solute_box_gro_path: str,
        input_top_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        additional_notes: Optional[str] = None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS solvate",
            "details": f"added solute at {solute_box_gro_path} to {solvent_box_gro_path}",
            "action(s)": f"saved at {output_base_dir}/{run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.OUTPUT_GRO_NAME}, saved raw topology file at {output_base_dir}/{run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.TOPOL_NAME}",
            "additional_notes": additional_notes,
        }
