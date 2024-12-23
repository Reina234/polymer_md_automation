from A_modules.shared.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory

from A_modules.shared.packmol.solvent_box import PackmolSolventBox
from data_models.solvent import Solvent
import os
from A_modules.shared.file_conversion.converters.editconf_gro_to_pdb import (
    EditconfGROtoPDBConverter,
)
from A_modules.shared.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)
from A_modules.atomistic.gromacs.equilibriation.mdp_cache import MDPCache
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from A_modules.atomistic.gromacs.equilibriation.workflow_step.npt_workflow_step import (
    NptWorkflowStep,
)
from A_modules.atomistic.gromacs.commands.grompp import Grompp
from A_modules.atomistic.gromacs.commands.mdrun import MDrun

from A_modules.shared.utils.file_utils import check_directory_exists


test = ConverterFactory().get_converter("pdb", "mol2")
mol2_file = test.run("input/solvents/pdb/hexane.pdb", "TEST", verbose=True)
import pandas as pd

print("!!!" + mol2_file)
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)

parameterizer = ACPYPEParameterizer()
file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)

parametized = parameterizer.run(mol2_file, "TEST", file_config, verbose=True)

from A_modules.shared.pdb_validation.gromacs_validator import GROMACSPDBValidator

updated_path = "input/solvents/pdb/hexane.pdb"

# gromacsvalidator = GROMACSPDBValidator()
# updated_path = "TEST/edited_pdb.pdb"
# gromacsvalidator.validate(
#    "input/solvents/pdb/hexane.pdb",
#    output_file_path=updated_path,
#    padding=0.5 * 10**-10,
# )
gro = parametized.gro_path
top = parametized.top_path


converter = EditconfGROtoPDBConverter()
box_size_nm = [5, 5, 5]
# may need to fix file extensions on converters
output_pdb = converter.run(gro, "ZZZZZ", output_name="test_output.pdb")


# Path to a valid solvent PDB file
output_dir = "test_output"
output_file = "test3_solvent_box.pdb"
# Example box size in nanometers
solvent = Solvent("Hexane", 86, 660, output_pdb, "TMZK")

# Create an instance of PackmolSolventBox
packmol_operation = PackmolSolventBox()

# Run the Packmol operation
result = packmol_operation.run(
    output_pdb,
    output_file,
    output_dir="ZZZZZ",
    solvent=solvent,
    box_size_nm=box_size_nm,
)

# Validate the result
assert os.path.exists(result), f"Output file not found: {result}"

converter = EditconfPDBtoGROConverter()
output_gro = converter.run(
    result, "ZZZZZ", box_size_nm=box_size_nm, output_name="test_output.gro"
)

print(result)

# Initialize MDP cache
mdp_cache = MDPCache(cache_dir="cache")

# Set up workflow
workflow = FullEquilibrationWorkflow(mdp_cache)

workflow.add_step(
    workflow_step=NptWorkflowStep(Grompp(), MDrun()),
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/em.mdp",
    base_params={"nsteps": "50000"},
)

workflow.add_step(
    workflow_step=NptWorkflowStep(Grompp(), MDrun()),
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/em_2.mdp",
    base_params={"nsteps": "50000"},
)


workflow.add_step(
    workflow_step=NptWorkflowStep(Grompp(), MDrun()),
    template_path="A_modules/atomistic/gromacs/equilibriation/templates/npt.mdp",
    base_params={"pressure": "1.0", "compressibility": "4.5e-5"},
)

# List of varying parameters
varying_params_list = [{"temp": str(temp), "pressure_tau": "2.0"} for temp in [300]]

# Run workflow
workflow.run(
    input_gro_path=output_gro,
    input_topol_path="temp/TMZK.acpype/TMZK_GMX.top",
    temp_output_dir="temp_outputs",
    main_output_dir="final_outputs",
    keep_files=["gro", "log"],
    varying_params_list=varying_params_list,
    verbose=True,
)
