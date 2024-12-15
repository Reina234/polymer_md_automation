from A_modules.shared.file_conversion.converters.obabel_pdb_to_mol2_converter import (
    OBabelPDBtoMOL2Converter,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory

converter = OBabelPDBtoMOL2Converter()
converter.run("styrene.pdb", "TEST", verbose=True)

from A_modules.shared.utils.utils import check_directory_exists

test = ConverterFactory().get_converter("pdb", "mol2")
mol2_file = test.run("styrene.pdb", "TEST", verbose=True)

print("!!!" + mol2_file)
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)

parameterizer = ACPYPEParameterizer()
file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)

parameterizer.run(mol2_file, "TEST", file_config, verbose=True)
