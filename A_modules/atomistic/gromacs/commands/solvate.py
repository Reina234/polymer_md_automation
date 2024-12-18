from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.file_utils import (
    file_exists_check_wrapper,
    directory_exists_check_wrapper,
)
from typing import Optional, List, Dict
import os


class Solvate(BaseGromacsCommand):

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    @file_exists_check_wrapper(file_arg_index=1)
    @file_exists_check_wrapper(file_arg_index=2)
    @file_exists_check_wrapper(file_arg_index=3)
    @directory_exists_check_wrapper(dir_arg_index=4)
    def run(
        self,
        solute_gro_path: str,
        solvent_gro_path: str,
        input_topol_path: str,
        output_dir: str,
        output_name: str = "solvated.gro",
        additional_notes: Optional[str] = None,
    ):
        output_gro_path = os.path.join(output_dir, output_name)

        command = self._create_solvate_command(
            solute_gro_path=solute_gro_path,
            solvent_gro_path=solvent_gro_path,
            input_topol_path=input_topol_path,
            output_gro_path=output_gro_path,
        )

        self._execute(command)
        if self.metadata_tracker:
            self.metadata_tracker.update_metadata(
                self.metadata(
                    solute_gro_path=solute_gro_path,
                    solvent_gro_path=solvent_gro_path,
                    input_topol_path=input_topol_path,
                    output_gro_path=output_gro_path,
                    additional_notes=additional_notes,
                )
            )
        return output_gro_path

    def _create_solvate_command(
        self,
        solute_gro_path: str,
        solvent_gro_path: str,
        input_topol_path: str,
        output_gro_path: str,
    ) -> List[str]:

        solvate_command = [
            "gmx",
            "solvate",
            "-cp",
            solute_gro_path,
            "-cs",
            solvent_gro_path,
            "-p",
            input_topol_path,
            "-o",
            output_gro_path,
        ]
        return solvate_command

    def metadata(
        self,
        solute_gro_path: str,
        solvent_gro_path: str,
        input_topol_path: str,
        output_gro_path: str,
        additional_notes: Optional[str],
    ) -> Dict[str, str]:
        return {
            "program(s) used": "GROMACS solvate",
            "details": f"added solute at {solute_gro_path} to {solvent_gro_path}",
            "action(s)": f"saved at {output_gro_path}, topol at {input_topol_path}",
            "additional_notes": additional_notes,
        }
