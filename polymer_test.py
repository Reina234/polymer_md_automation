from data_models.solvent import Solvent
from A_modules.atomistic.gromacs.utils.utils import create_solvent_box_gro
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory
from A_modules.atomistic.gromacs.workflows.equilibriated_solvent_box.file_preparation_utils import (
    process_solvent_files,
)
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from typing import List, Optional
from data_models.output_types import GromacsPaths
from A_config.paths import TEMP_DIR, LOG_DIR
from A_modules.shared.utils.file_utils import add_identifier_name
from A_modules.atomistic.gromacs.utils.mdp_utils import format_temperatures
from A_modules.atomistic.gromacs.commands.insert_molecules import InsertMolecules
import os

# NOTE: need to add cleanup for temp files

temp_dir_run = "temp_29_test"
polymer_pdb = "styrene.pdb"
parameterizer = ACPYPEParameterizer()
mol2_converter = ConverterFactory().get_converter("pdb", "mol2")
mol2_file = mol2_converter.run(polymer_pdb, temp_dir_run)
file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
parametized_files = parameterizer.run(mol2_file, temp_dir_run, file_config)

solvent_equilibriated_gro = "TEST_RUN_28/hexane_run_1/equilibriated_gros/temp_298.gro"
insert_molecules = InsertMolecules()
# run insert-molecules
solute_in_box = insert_molecules.run(
    solvent_equilibriated_gro, parametized_files.gro_path, 1, temp_dir_run, "solute_box"
)
