from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_modules.shared.utils.utils import file_exists_check_wrapper
from A_config.constants import MassUnits2

class InsertMolecules(BaseGromacsCommand):
    default_mass_units = MassUnits2.GRAM

    def __init__(self, metadata_tracker: MetadataTracker = None):
        super().__init__(metadata_tracker)

    
    def run(
            self, 
            input_pdb_path: str,
            gro_path: str,
    )
