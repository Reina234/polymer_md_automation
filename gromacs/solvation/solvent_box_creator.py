import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import GROMACS_OUTPUT_SUBDIR, BASE_OUTPUT_DIR
from config.constants import DEFAULT_BOX_SIZE_NM
from typing import Optional, List
from preprocessing.calculation_utils import calculate_num_particles
from config.constants import LengthUnits
from gromacs.base_gromacs_command import BaseGromacsCommand

# NOTE: add in file type checks


class SolventInsertion(BaseGromacsCommand):
    OUTPUT_NAME = "solvent_box.gro"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def run(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        additional_notes: Optional[str] = None,
    ) -> str:
        command, output_gro_path = self._create_editconf_command(
            solvent_pdb_path=solvent_pdb_path,
            run_name=run_name,
            solvent_density=solvent_density,
            solvent_molecular_weight=solvent_molecular_weight,
            output_base_dir=output_base_dir,
            box_size_nm=box_size_nm,
            additional_notes=additional_notes,
        )
        self._execute(command)
        if self.metadata_tracker:
            self._update_metadata(
                solvent_pdb_path=solvent_pdb_path,
                run_name=run_name,
                solvent_density=solvent_density,
                solvent_molecular_weight=solvent_molecular_weight,
                output_base_dir=output_base_dir,
                box_size_nm=box_size_nm,
                additional_notes=additional_notes,
            )
        return output_gro_path

    def _create_editconf_command(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        additional_notes: Optional[str] = None,
    ) -> str:
        num_particles = calculate_num_particles(
            box_size_nm,
            solvent_molecular_weight,
            solvent_density,
            box_units=self.UNITS,
        )
        solvent_box_gro_path = os.path.join(
            output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR, self.OUTPUT_NAME
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
            solvent_box_gro_path,
        ]
        return editconf_command, solvent_box_gro_path

    def metadata(
        self,
        solvent_pdb_path: str,
        run_name: str,
        solvent_density: float,
        solvent_molecular_weight: float,
        output_base_dir: str = BASE_OUTPUT_DIR,
        box_size_nm: List[float] = DEFAULT_BOX_SIZE_NM,
        additional_notes: Optional[str] = None,
    ) -> dict:
        return {
            "program(s) used": "GROMACS insert-molecules",
            "details": f"created a solvent box of size {box_size_nm} with units {self.UNITS.value}",
            "action(s)": f"used molecule at {solvent_pdb_path}, saved at {output_base_dir}/{run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.OUTPUT_NAME}",
            "additional_notes": additional_notes,
        }
