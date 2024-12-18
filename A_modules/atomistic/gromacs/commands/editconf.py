from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.file_utils import (
    file_exists_check_wrapper,
    directory_exists_check_wrapper,
)
from typing import Optional, List, Dict
import os


class Editconf(BaseGromacsCommand):
    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    @file_exists_check_wrapper(file_arg_index=1)
    @directory_exists_check_wrapper(dir_arg_index=2)
    def run(
        self,
        input_gro_or_pdb_path: str,
        output_dir: str,
        box_size_nm: Optional[List[float]] = None,
        output_name: str = "edited_box.gro",
        additional_notes: Optional[str] = None,
    ):
        output_gro_path = os.path.join(output_dir, output_name)

        command = self._create_editconf_command(
            input_gro_or_pdb_path=input_gro_or_pdb_path,
            output_gro_path=output_gro_path,
            box_size_nm=box_size_nm,
        )

        self._execute(command)
        if self.metadata_tracker:
            self.metadata_tracker.update_metadata(
                self.metadata(
                    input_gro_path=input_gro_or_pdb_path,
                    output_gro_path=output_gro_path,
                    box_size_nm=box_size_nm,
                    additional_notes=additional_notes,
                )
            )
        return output_gro_path

    def _create_editconf_command(
        self,
        input_gro_or_pdb_path: str,
        output_gro_path: str,
        box_size_nm: Optional[List[float]],
    ) -> List[str]:
        command = [
            "gmx",
            "editconf",
            "-f",
            input_gro_or_pdb_path,
            "-o",
            output_gro_path,
        ]
        if box_size_nm is not None:
            command.extend(
                [
                    "-box",
                    str(box_size_nm[0]),
                    str(box_size_nm[1]),
                    str(box_size_nm[2]),
                ]
            )
        return command

    def metadata(
        self,
        input_gro_path: str,
        output_gro_path: str,
        box_size_nm: Optional[List[float]],
        additional_notes: Optional[str],
    ) -> Dict[str, str]:
        details = f"created a polymer box"
        if box_size_nm:
            details += f" of size {box_size_nm} with units {self.default_units}"
        return {
            "program(s) used": "GROMACS - editconf",
            "details": details,
            "action(s)": f"used molecule at {input_gro_path}, saved at {output_gro_path}",
            "additional_notes": additional_notes,
        }
