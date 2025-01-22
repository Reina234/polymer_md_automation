from rdkit_new.homopolymer_generator import HomopolymerGenerator
from config.paths import TEMP_DIR
from config.acpype_config import AcpypeOutputConfig
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from modules.shared.file_conversion.converter_factory import ConverterFactory
from modules.shared.utils.file_utils import copy_file, delete_directory
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from data_models.output_types import GromacsPaths
from config.paths import PARAMETERISED_POLYMER_DIR
import os

# NOTE: I think I can use just a copolymer generator, and create a homopolymer by passing in two of the same monomer units


def run_polymerisation_workflow(
    monomer_smiles: str,
    num_units: int,
    output_dir: str = PARAMETERISED_POLYMER_DIR,
    temp_dir: str = TEMP_DIR,
    keep_pdb: bool = False,
    keep_mol2: bool = False,
    uff_optimise: bool = True,
    overwrite: bool = True,
    polymer_name: str = "POLY",
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    converter_factory: ConverterFactory = ConverterFactory(),
    cleanup: bool = True,
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
    if cleanup:
        delete_directory(temp_dir, confirm=False)
    return parametized_polymer
