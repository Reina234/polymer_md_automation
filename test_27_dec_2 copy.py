from data_models.solvent import Solvent
from A_modules.atomistic.gromacs.utils.utils import create_solvent_box_gro
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory
from A_modules.atomistic.gromacs.utils.utils import (
    process_solvent_files,
)


test_dir = "TEST_RUN_27_2"
acpype_subdir = "acpype_output"


solvent = Solvent("Hexane", 86, 660, "input/solvents/pdb/hexane.pdb", "TMZK")
box_size_nm = [5, 5, 5]


mol2_converter = ConverterFactory().get_converter("pdb", "mol2")
mol2_file = mol2_converter.run(solvent.pdb_path, test_dir, verbose=True)
parameterizer = ACPYPEParameterizer()
file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
parametized = parameterizer.run(mol2_file, test_dir, file_config, verbose=True)
gro = parametized.gro_path
top = parametized.top_path
itp = parametized.itp_path


solvent_box_gro = create_solvent_box_gro(gro, test_dir, box_size_nm, solvent)
output_files = process_solvent_files(
    top, itp, solvent_box_gro, solvent.pdb_molecule_name, test_dir, test_dir
)
