import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from typing import Optional, List
from preprocessing.pdb_utils import calculate_num_particles
from config.constants import LengthUnits

# NOTE: add in file type checks


class SolventPlacer:
    OUTPUT_GRO_NAME = "solvent_box.gro"
    UNITS = LengthUnits.NANOMETER

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        self.metadata_tracker = metadata_tracker

    def create_box(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = [3.0, 3.0, 3.0],
        additional_notes: Optional[str] = None,
    ) -> str:
        num_particles = calculate_num_particles(
            box_size_nm,
            solvent_molecular_weight,
            solvent_density,
            box_units=self.UNITS,
        )
        output_file_path = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, self.OUTPUT_GRO_NAME
        )

        output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
        os.makedirs(output_dir, exist_ok=True)
        editconf_command = [
            "gmx",
            "insert-molecules",
            "-ci",
            solvent_pdb_path,
            "-nmol",
            str(round(num_particles)),
            "-box",
            str(box_size_nm[0]),
            str(box_size_nm[1]),
            str(box_size_nm[2]),
            "-o",
            output_file_path,
        ]
        subprocess.run(editconf_command, check=True)
        if self.metadata_tracker:
            self._update_metadata(
                solvent_pdb_path, output_file_path, box_size_nm, additional_notes
            )
        return output_file_path

    def metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        box_size: List[float],
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
        box_size: List[float],
        additional_notes: Optional[str] = None,
    ) -> None:
        metadata = self.metadata(
            input_file_path, output_file_path, box_size, additional_notes
        )
        self.metadata_tracker.add_step(step_name="GROMACS", details=metadata)
