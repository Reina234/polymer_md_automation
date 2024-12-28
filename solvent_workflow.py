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
from A_modules.atomistic.gromacs.equilibriation.full_equilibriation_workflow import (
    FullEquilibrationWorkflow,
)
from typing import List, Dict, Optional
from A_modules.atomistic.config import GromacsPaths
from A_config.paths import TEMP_DIR, LOG_DIR


def add_identifier_name(
    name: str,
    default_identifier: Optional[str] = None,
    identifier: Optional[str] = None,
    suffix: Optional[str] = None,
) -> str:
    parts = [name.lower()]
    if default_identifier:
        parts.append(default_identifier)
    if identifier:
        parts.append(identifier)
    if suffix:
        parts.append(suffix)
    return "_".join(parts)


def format_temperatures(temperatures: List[float]) -> List[Dict[str, str]]:
    return [{"temp": str(temp)} for temp in temperatures]


def solvent_workflow(
    solvent: Solvent,
    box_size_nm: List[float],
    output_dir: str,
    workflow: FullEquilibrationWorkflow,
    temperatures: List[float],
    identifier: Optional[str] = None,
    temp_dir=TEMP_DIR,
    log_dir=LOG_DIR,
    verbose: bool = True,
    save_intermediate_edr=True,
    save_intermediate_gro=True,
    save_intermediate_log=True,
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
) -> GromacsPaths:
    parameterizer = parameterizer(solvent.pdb_molecule_name)
    itp_name = add_identifier_name(solvent.name, identifier=identifier, suffix=None)
    mol2_converter = ConverterFactory().get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(solvent.pdb_path, temp_dir, verbose=verbose)
    file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
    parametized_files = parameterizer.run(
        mol2_file, temp_dir, file_config, verbose=verbose
    )
    solvent_box_gro = create_solvent_box_gro(
        parametized_files.gro_path, temp_dir, box_size_nm=box_size_nm, solvent=solvent
    )
    reformatted_files = process_solvent_files(
        parametized_files.itp_path,
        solvent_box_gro,
        parametized_files.top_path,
        new_residue_name=solvent.pdb_molecule_name,
        output_itp_dir=output_dir,
        output_gro_dir=temp_dir,
        output_topol_dir=temp_dir,
        output_itp_name=itp_name,
    )

    temperatures = format_temperatures(temperatures)

    gro_name = add_identifier_name(
        solvent.name,
        default_identifier="equilibriated",
        identifier=identifier,
        suffix=None,
    )
    equilibriated_gro = workflow.run(
        reformatted_files.gro_path,
        reformatted_files.top_path,
        temp_dir,
        output_dir,
        log_dir,
        final_gro_name=gro_name,
        varying_params_list=temperatures,
        save_intermediate_edr=save_intermediate_edr,
        save_intermediate_gro=save_intermediate_gro,
        save_intermediate_log=save_intermediate_log,
        verbose=verbose,
    )

    return GromacsPaths(
        itp_path=reformatted_files.itp_path, top_path=None, gro_path=equilibriated_gro
    )
