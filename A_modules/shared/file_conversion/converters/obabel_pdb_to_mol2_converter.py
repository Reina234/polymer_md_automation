import os
import subprocess
from A_modules.shared.file_conversion.converters.base_converter import (
    BaseConverter,
)
from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.shared.utils.utils import construct_output_file_path
from typing import Optional, List, Tuple


class OBabelPDBtoMOL2Converter(BaseConverter):
    input_file_type = "pdb"
    output_file_type = "mol2"
    program = "Open Babel"

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    def _create_obabel_command(
        self, input_file_path: str, output_dir: Optional[str] = None
    ) -> Tuple[List, str]:
        output_file_path = construct_output_file_path(
            file_path=input_file_path,
            new_output_dir=output_dir,
            new_file_extension=self.output_file_type,
            suppress_warning=True,
            expected_file_type=self.input_file_type,
        )
        command = [
            "obabel",
            input_file_path,
            "-O",
            output_file_path,
            "--gen3D",
            "--h",  # add 3D coords via --h
        ]
        return command, output_file_path

    def run(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ):
        command, output_file_path = self._create_obabel_command(
            input_file_path, output_dir
        )
        self._execute(command)

        if self.metadata_tracker is not None:
            self._update_metadata(input_file_path, output_file_path, additional_notes)

        return output_file_path
