from abc import ABC, abstractmethod
from typing import Dict, List, Optional
import os
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun
from A_modules.atomistic.gromacs.utils.mdp_utils import generate_mdp_file
import logging
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.file_utils import cleanup_directory, create_temp_directory
from A_config.paths import TEMP_DIR
from A_modules.shared.metadata_tracker import MetadataTracker

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


# NOTE: pass in mdp later, consider using a dict to decide which params to replace, maybe use dependecy injection to pass in what template to use at the workflow level


class BaseEquilibrationWorkflow(ABC):
    default_name: str

    def __init__(
        self,
        grompp: BaseGromacsCommand = Grompp,
        mdrun: BaseGromacsCommand = MDrun,
        metadata_tracker: MetadataTracker = None,
    ):
        self.grompp = grompp(metadata_tracker)
        self.mdrun = mdrun(metadata_tracker)

    @abstractmethod
    @property
    def mdp(self) -> str:
        pass

    @abstractmethod
    @mdp.setter
    def mdp(self, **kwargs: str) -> str:
        pass

    def run(
        self,
        mdp_file_path: str,
        input_gro_path: str,
        input_topol_path: str,
        output_dir: str = TEMP_DIR + "/equilibration",
        output_name: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ) -> str:
        if output_name is None:
            output_name = self.default_name

        grompp_output = self.grompp.run(
            mdp_file_path,
            input_gro_path,
            input_topol_path,
            output_tpr_dir=output_dir,
            output_tpr_name=output_name,
            additional_notes=additional_notes,
            verbose=verbose,
        )
        mdrun_output = self.mdrun.run(
            grompp_output,
            output_name,
            additional_notes=additional_notes,
            verbose=verbose,
        )
        return mdrun_output
