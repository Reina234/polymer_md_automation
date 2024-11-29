import os
import subprocess
import shutil
import logging
from typing import Optional, Dict, List
from preprocessing.parameterizers.base_parameterizer import BaseParameterizer
from preprocessing.metadata_tracker import MetadataTracker
from preprocessing.utils import (
    check_file_exists,
    check_file_type,
    copy_files,
    rename_basenames,
)
from config.paths import (
    ACPYPE_PARAMETERIZER_OUTPUT_SUBDIR,
    ACPYPE_BASE_NAME,
    TEMPORARY_OUTPUT_DIR,
    BASE_OUTPUT_DIR,
    ACPYPE_SOLVENT_OUTPUT_SUBDIR,
    FilesToExtract,
)
from enum import Enum


logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class ACPYPEParameterizer(BaseParameterizer):
    """
    Parameterizes a .mol2 polymer using ACPYPE for atomistic simulations.
    Generates directories containing .gro, .top, and .itp files.
    """

    def __init__(
        self,
        metadata_tracker: Optional[MetadataTracker] = None,
        acpype_molecule_name: str = ACPYPE_BASE_NAME,
    ):
        super().__init__(metadata_tracker)
        self.molecule_name = acpype_molecule_name
        self.raw_output_dir = os.path.join(
            TEMPORARY_OUTPUT_DIR, f"{self.molecule_name}.acpype"
        )

    def parameterize(
        self,
        input_file_path: str,
        output_dir: str,
        files_to_extract: List[FilesToExtract],
        rename_output: Optional[str] = None,
        posre: bool = False,
        additional_notes: Optional[str] = None,
    ) -> str:

        abs_input_file_path, abs_output_dir = self._prepare_directories(
            input_file_path, output_dir
        )
        self._run_acpype_command(abs_input_file_path, TEMPORARY_OUTPUT_DIR)

        files_to_move = [
            f"{self.molecule_name}_GMX.{file_type.value}"
            for file_type in files_to_extract
        ]

        if posre:
            files_to_move.append(f"posre_{self.molecule_name}.itp")

        moved_files = copy_files(files_to_move, self.raw_output_dir, abs_output_dir)

        if rename_output:
            moved_files = rename_basenames(moved_files, rename_output)

        if self.metadata_tracker is not None:
            self._update_metadata(abs_input_file_path, abs_output_dir, additional_notes)

        logger.info("ACPYPE parameterization completed successfully.")
        return moved_files

    def _prepare_directories(self, input_file_path: str, output_dir: str) -> None:
        check_file_exists(input_file_path)
        check_file_type(input_file_path, "mol2")
        os.makedirs(output_dir, exist_ok=True)
        abs_input_file_path_full = os.path.abspath(input_file_path)
        abs_output_dir_full = os.path.abspath(output_dir)
        return abs_input_file_path_full, abs_output_dir_full

    def _run_acpype_command(self, input_file_path: str, temp_dir: str):
        """
        Runs the ACPYPE command in the temporary directory.
        """
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

        logger.info("Running ACPYPE...")
        logger.debug(f"ACPYPE command: {' '.join(acpype_command)}")

        # Overwrite .acpype directory if it exists
        acpype_temp_dir = os.path.join(temp_dir, f"{self.molecule_name}.acpype")
        if os.path.exists(acpype_temp_dir):
            shutil.rmtree(acpype_temp_dir)
        os.makedirs(acpype_temp_dir, exist_ok=True)

        result = subprocess.run(
            acpype_command,
            cwd=temp_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        logger.debug(f"ACPYPE stdout: {result.stdout}")
        logger.debug(f"ACPYPE stderr: {result.stderr}")

        if result.returncode != 0:
            raise RuntimeError(
                f"ACPYPE failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}"
            )

    def metadata(
        self,
        input_file_path: str,
        output_dir: str,
        additional_notes: Optional[str] = None,
    ) -> Dict:
        """
        Returns metadata about ACPYPE.
        """
        return {
            "program(s) used": "ACPYPE (GAFF2)",
            "details": "generates GROMACS-compatible force field files. Using GAFF2",
            "action(s)": f"parameterized {input_file_path}, saved to {output_dir}",
            "additional_note(s)": additional_notes,
        }
