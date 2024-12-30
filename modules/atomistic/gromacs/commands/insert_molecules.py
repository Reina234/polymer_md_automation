from modules.shared.metadata_tracker import MetadataTracker
from modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    prepare_output_file_path,
)
from config.constants import MassUnits
from typing import Optional, Dict, Tuple, List


class InsertMolecules(BaseGromacsCommand):
    default_mass_units = MassUnits.GRAM

    def __init__(self, metadata_tracker: MetadataTracker = None):
        super().__init__(metadata_tracker)

    @file_type_check_wrapper(file_arg_index=1, expected_file_type="gro")
    @file_type_check_wrapper(file_arg_index=2, expected_file_type="gro")
    def run(
        self,
        input_box_gro_path: str,
        solvent_gro_path: str,
        num_molecules: int,
        output_dir: Optional[str] = None,
        output_name: Optional[str] = None,
        additional_notes: Optional[str] = None,
    ):
        output_path = prepare_output_file_path(
            input_box_gro_path,
            output_extension="gro",
            output_dir=output_dir,
            output_name=output_name,
        )
        command, output_path = self._create_command(
            input_box_gro_path=input_box_gro_path,
            solvent_gro_path=solvent_gro_path,
            num_molecules=num_molecules,
            output_path=output_path,
        )

        self._execute(command)
        if self.metadata_tracker:
            self._update_metadata(
                input_box_gro_path=input_box_gro_path,
                solvent_gro_path=solvent_gro_path,
                num_molecules=num_molecules,
                output_path=output_path,
                additional_notes=additional_notes,
            )

        return output_path

    def _create_command(
        self,
        input_box_gro_path: str,
        solvent_gro_path: str,
        num_molecules: int,
        output_path: str = None,
    ) -> Tuple[List[str], str]:

        if not output_path:
            output_path = input_box_gro_path
            append_flag = "-append"
        else:
            append_flag = None

        # Build the GROMACS command
        command = [
            "gmx",
            "insert-molecules",
            "-f",
            input_box_gro_path,  # Input structure file with the existing solvent
            "-ci",
            solvent_gro_path,  # Input solvent file to insert
            "-nmol",
            str(num_molecules),  # Number of molecules to insert
            "-o",
            output_path,  # Output file
        ]

        # Add the append flag if applicable
        if append_flag:
            command.append(append_flag)

        return command, output_path

    def metadata(
        self,
        input_box_gro_path: str,
        solvent_gro_path: str,
        num_molecules: int,
        output_path: str,
        additional_notes: Optional[str],
    ) -> Dict[str, str]:
        return {
            "program(s) used": "GROMACS",
            "command": "insert-molecules",
            "details": f"Insert {num_molecules} molecules from {solvent_gro_path} into {input_box_gro_path}",
            "additional_notes": additional_notes,
        }
