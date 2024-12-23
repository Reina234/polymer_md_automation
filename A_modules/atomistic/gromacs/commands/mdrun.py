from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    prepare_output_file_path,
)
from A_config.constants import MassUnits2
from typing import Optional, Dict, Tuple, List
import os


# NOTE: NEED TO FIX, think abt temp, and where to input/output dir - e.g. in the new base class for equilibrium or not
class MDrun(BaseGromacsCommand):
    def __init__(self, metadata_tracker=None):
        super().__init__(metadata_tracker)

    def run(
        self,
        input_tpr_path: str,
        output_name: str,
        verbose: bool = False,
        additional_flags: Optional[List[str]] = None,
        additional_notes: Optional[str] = None,
    ) -> Dict[str, str]:
        command = self._create_command(input_tpr_path, output_name, additional_flags)
        self._execute(command, verbose=verbose)
        output_files = {
            ext: f"{output_name}.{ext}" for ext in ["gro", "log", "edr", "trr"]
        }
        if self.metadata_tracker:
            self.metadata_tracker.update_metadata(
                self.metadata(input_tpr_path, output_name, additional_notes)
            )
        return {k: v for k, v in output_files.items() if os.path.exists(v)}

    def _create_command(
        self,
        input_tpr_path: str,
        output_name: str,
        additional_flags: Optional[List[str]] = None,
    ) -> List[str]:
        command = ["gmx", "mdrun", "-s", input_tpr_path, "-deffnm", output_name]
        if additional_flags:
            command.extend(additional_flags)
        return command

    def metadata(
        self, input_tpr_path: str, output_name: str, additional_notes: Optional[str]
    ):
        return {
            "input_tpr_path": input_tpr_path,
            "output_name": output_name,
            "additional_notes": additional_notes,
        }
