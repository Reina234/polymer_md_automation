from modules.shared.metadata_tracker import MetadataTracker
from modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    check_file_type,
    prepare_output_file_path,
)
from config.constants import MassUnits
from typing import Optional, Dict, Tuple, List


class GenIon(BaseGromacsCommand):
    default_mass_units = MassUnits.GRAM

    def __init__(self, metadata_tracker: MetadataTracker = None):
        super().__init__(metadata_tracker)

    def _create_command(
        self,
        input_box_gro_path: str,
        tpr_path: str,
        top_path: str,
        output_path: str,
        pname: str,
        nname: str,
        conc: Optional[float],
    ) -> Tuple[List[str], str]:

        command = [
            "gmx",
            "genion",
            "-s",
            tpr_path,  # Number of molecules to insert
            "-o",
            output_path,  # Output file
            "-p",
            top_path,
            "-pname",
            pname,
            "-nname",
            nname,
            "-neutral",
        ]

        # Add the append flag if applicable
        if conc:
            command.append("-conc")
            command.append(str(conc))

        return command, output_path

    def run(
        self,
        input_box_gro_path: str,
        tpr_path: str,
        top_path: str,
        output_dir: Optional[str] = None,
        output_name: str = "neutralised_polymer.gro",
        pname: str = "NA",
        nname: str = "CL",
        conc: Optional[float] = None,
        verbose: bool = True,
    ):
        check_file_type(input_box_gro_path, "gro")
        check_file_type(tpr_path, "tpr")
        check_file_type(top_path, "top")

        output_path = prepare_output_file_path(
            input_box_gro_path,
            output_extension="gro",
            output_dir=output_dir,
            output_name=output_name,
        )
        command, output_path = self._create_command(
            input_box_gro_path=input_box_gro_path,
            tpr_path=tpr_path,
            top_path=top_path,
            output_path=output_path,
            pname=pname,
            nname=nname,
            conc=conc,
        )

        self._execute(command, verbose=verbose)
        return output_path

    def metadata(self):
        return None
