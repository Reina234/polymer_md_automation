from A_modules.shared.file_conversion.converter_factory import ConverterFactory
from A_modules.shared.utils.file_utils import (
    overwrite_directory,
    file_exists_check_wrapper,
    directory_exists_check_wrapper,
)
from A_modules.shared.command_line_operation import CommandLineOperation
from A_modules.atomistic.acpype_parameterizer.acpype_config import (
    AcpypeOutputConfig,
)
from data_models.output_types import GromacsPaths
from A_modules.atomistic.acpype_parameterizer.acpype_utils import (
    generate_acpype_paths,
    copy_acpype_files,
    rename_acpype_paths,
)
from A_config.paths import ACPYPE_POLYMER_NAME, TEMP_DIR
import os
from typing import Optional


class ACPYPEParameterizer(CommandLineOperation):

    def __init__(
        self,
        metadata_tracker=None,
        acpype_molecule_name: str = ACPYPE_POLYMER_NAME,
        temp_dir: str = TEMP_DIR,
    ):
        self.molecule_name = acpype_molecule_name
        self._temp_dir = temp_dir

        super().__init__(metadata_tracker)

    @property
    def temp_dir(self):
        return self._temp_dir

    @temp_dir.setter
    def temp_dir(self, value):
        self._temp_dir = value
        # Automatically update dependent attributes
        self.raw_output_dir = os.path.join(value, f"{self.molecule_name}.acpype")

    @property
    def raw_output_dir(self):
        return os.path.join(self.temp_dir, f"{self.molecule_name}.acpype")

    def acpype_command(self, input_file_path: str):
        acpype_command = [
            "acpype",
            "-i",
            input_file_path,
            "-o",
            "gmx",
            "-n",
            "0",
            "-a",
            "gaff2",
            "-b",
            self.molecule_name,
        ]

        return acpype_command

    # NOTE: may need to allow for non mol_2 types
    @file_exists_check_wrapper(file_arg_index=1)
    @directory_exists_check_wrapper(dir_arg_index=2)
    def run(
        self,
        input_file_path: str,
        output_dir: str,
        acpype_output_config: AcpypeOutputConfig,
        new_file_name: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ) -> GromacsPaths:
        input_file_path = os.path.abspath(input_file_path)
        command = self.acpype_command(input_file_path=input_file_path)
        overwrite_directory(self.raw_output_dir)
        self._execute(command, cwd=self.temp_dir, verbose=verbose)

        raw_acpype_paths = generate_acpype_paths(
            acpype_output_config, self.raw_output_dir, self.molecule_name
        )
        copied_acpype_paths = copy_acpype_files(raw_acpype_paths, output_dir)
        if new_file_name:
            renamed_acpype_paths = rename_acpype_paths(
                copied_acpype_paths, new_file_name
            )

            final_acpype_paths = renamed_acpype_paths
        else:
            final_acpype_paths = copied_acpype_paths

        if self.metadata_tracker:
            self._update_metadata(input_file_path, final_acpype_paths, additional_notes)

        return final_acpype_paths

    def metadata(
        self,
        input_file_path,
        output_acpype_paths: GromacsPaths,
        additional_notes: Optional[str],
    ):
        return {
            "program(s) used": "ACPYPE",
            "details": f"parameterized {input_file_path} with ACPYPE",
            "action(s)": f"saved output files at {output_acpype_paths.to_list()}",
        }
