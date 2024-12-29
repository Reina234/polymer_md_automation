from A_modules.atomistic.rdkit.homopolymer_generator import HomopolymerGenerator
from A_config.paths import TEMP_DIR
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory
from A_modules.shared.utils.file_utils import copy_file

# NOTE: I think I can use just a copolymer generator, and create a homopolymer by passing in two of the same monomer units


def run_polymerisation_workflow(
    monomer_smiles: str,
    num_units: int,
    output_dir: str,
    temp_dir: str = TEMP_DIR,
    keep_pdb: bool = False,
    uff_optimise: bool = True,
    overwrite: bool = True,
    converter_factory: ConverterFactory = ConverterFactory(),
) -> str:
    generator = HomopolymerGenerator()
    polymer_pdb = generator.generate_polymer(
        monomer_smiles,
        num_units,
        output_dir=temp_dir,
        uff_optimise=uff_optimise,
        overwrite=overwrite,
    )
    if keep_pdb:
        polymer_pdb = copy_file(generator.pdb_path, output_dir)

    mol2_converter = converter_factory.get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(polymer_pdb, output_dir)
    return mol2_file
