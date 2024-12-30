from A_modules.atomistic.rdkit.homopolymer_generator import HomopolymerGenerator
from A_config.paths import TEMP_DIR
from A_modules.atomistic.acpype_parameterizer.acpype_config import AcpypeOutputConfig
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from A_modules.shared.file_conversion.converter_factory import ConverterFactory
from A_modules.shared.utils.file_utils import copy_file
from A_modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from data_models.output_types import GromacsPaths
import os

# NOTE: I think I can use just a copolymer generator, and create a homopolymer by passing in two of the same monomer units


def run_polymerisation_workflow(
    monomer_smiles: str,
    num_units: int,
    output_dir: str,
    temp_dir: str = TEMP_DIR,
    keep_pdb: bool = False,
    keep_mol2: bool = False,
    uff_optimise: bool = True,
    overwrite: bool = True,
    polymer_name: str = "POLY",
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    converter_factory: ConverterFactory = ConverterFactory(),
) -> GromacsPaths:
    subdir = f"{monomer_smiles.lower()}_{num_units}"
    main_output_dir = os.path.join(output_dir, subdir)
    generator = HomopolymerGenerator()
    polymer_pdb = generator.generate_polymer(
        monomer_smiles,
        num_units,
        output_dir=temp_dir,
        uff_optimise=uff_optimise,
        overwrite=overwrite,
    )
    if keep_pdb:
        polymer_pdb = copy_file(generator.pdb_path, main_output_dir)

    if keep_mol2:
        mol2_output_dir = main_output_dir
    else:
        mol2_output_dir = temp_dir
    mol2_converter = converter_factory.get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(polymer_pdb, mol2_output_dir)
    parameterizer = parameterizer(acpype_molecule_name=polymer_name)
    file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
    parametized_polymer = parameterizer.run(mol2_file, main_output_dir, file_config)
    return parametized_polymer
