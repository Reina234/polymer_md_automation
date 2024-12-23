from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    prepare_output_file_path,
)
from A_config.constants import MassUnits2
from typing import Optional, Dict, Tuple, List
import os


class Grompp(BaseGromacsCommand):

    def __init__(self, metadata_tracker=None):
        super().__init__(metadata_tracker)

    def run(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ) -> str:
        output_tpr_path = prepare_output_file_path(
            input_gro_path,
            output_extension="tpr",
            output_dir=output_dir,
            output_name=output_name,
        )
        command = self._create_command(
            mdp_file_path,
            input_gro_path,
            input_topol_path,
            output_tpr_path,
        )
        self._execute(command, verbose=verbose)
        if self.metadata_tracker:
            self._update_metadata(
                mdp_file_path=mdp_file_path,
                input_gro_path=input_gro_path,
                input_topol_path=input_topol_path,
                output_tpr_path=output_tpr_path,
                additional_notes=additional_notes,
            )
        return output_tpr_path

    @file_type_check_wrapper(file_arg_index=1, expected_file_type="mdp")
    @file_type_check_wrapper(file_arg_index=2, expected_file_type="gro")
    @file_type_check_wrapper(file_arg_index=3, expected_file_type="top")
    @file_type_check_wrapper(file_arg_index=4, expected_file_type="tpr")
    def _create_command(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_tpr_path: str,
    ) -> List[str]:

        return [
            "gmx",
            "grompp",
            "-f",
            mdp_file_path,
            "-c",
            input_gro_path,
            "-p",
            input_topol_path,
            "-o",
            output_tpr_path,
        ]

    def metadata(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_tpr_path: str,
        additional_notes: Optional[str] = None,
    ) -> Dict:
        return {
            "program": "GROMACS",
            "mdp_file_path": mdp_file_path,
            "input_gro_path": input_gro_path,
            "input_topol_path": input_topol_path,
            "output_tpr_path": output_tpr_path,
            "additional_notes": additional_notes,
        }
