import os
import subprocess
from modules.shared.file_conversion.converters.base_converter import (
    BaseConverter,
)
from modules.shared.metadata_tracker import MetadataTracker
from modules.shared.utils.file_utils import (
    construct_output_file_path,
    file_type_check_wrapper,
)
from typing import Optional, List, Tuple
from modules.atomistic.gromacs.commands.editconf import Editconf
from pathlib import Path


class EditconfGROtoPDBConverter(BaseConverter):
    input_file_type = "gro"
    output_file_type = "pdb"
    program = "GROMACS"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def _run_impl(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
        box_size_nm: Optional[List[float]] = None,
        output_name: Optional[str] = None,
    ):
        if self.metadata_tracker:
            editconf = Editconf(self.metadata_tracker)
        else:
            editconf = Editconf()
        if output_name is None:
            output_name = Path(input_file_path).stem + ".pdb"
        output_path = editconf.run(
            input_file_path, output_dir, box_size_nm, output_name
        )
        return output_path
