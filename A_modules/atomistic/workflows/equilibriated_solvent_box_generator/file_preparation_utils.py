from A_modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    prepare_output_file_path,
)
from collections import OrderedDict
from A_modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from A_modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler
from typing import Optional
import logging
import os
from data_models.output_types import GromacsPaths
from A_modules.atomistic.gromacs.utils.utils import (
    get_gro_handler,
    get_residue_number,
    rename_residue_name_from_handler,
    rename_data_column_content,
    export_gro_handler,
    create_includes_section,
    delete_all_include_sections,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


@file_type_check_wrapper(file_arg_index=0, expected_file_type="itp")
@file_type_check_wrapper(file_arg_index=1, expected_file_type="gro")
@file_type_check_wrapper(file_arg_index=2, expected_file_type="top")
def process_solvent_files(
    input_itp_file: str,
    input_gro_file: str,
    input_top_file: str,
    forcefield: str = "amber99sb-ildn.ff/forcefield.itp",
    new_residue_name: Optional[str] = None,
    output_itp_dir: Optional[str] = None,
    output_gro_dir: Optional[str] = None,
    output_topol_dir: Optional[str] = None,
    output_itp_name: Optional[str] = None,
    output_gro_name: Optional[str] = None,
    output_topol_name: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
) -> GromacsPaths:
    if new_residue_name:
        if len(new_residue_name) > 5:
            raise ValueError(
                "Residue name must be 5 characters or less. Gromacs has fixed width for residue names."
            )
    output_itp_file = process_solvent_itp(
        input_itp_file=input_itp_file,
        new_residue_name=new_residue_name,
        output_dir=output_itp_dir,
        output_name=output_itp_name,
    )
    output_itp_path = os.path.abspath(output_itp_file)
    gro_handler = get_gro_handler(input_gro_file)
    residue_number = get_residue_number(gro_handler)
    output_top_file = prepare_solvent_topol(
        input_top_file=input_top_file,
        residue_number=residue_number,
        forcefield=forcefield,
        new_include_file=output_itp_path,
        new_residue_name=new_residue_name,
        output_name=output_topol_name,
        output_dir=output_topol_dir,
        parser=parser,
        del_posre=True,
        del_defaults=True,
    )
    if new_residue_name:
        gro_handler = rename_residue_name_from_handler(gro_handler, new_residue_name)

    output_gro_file = prepare_output_file_path(
        input_gro_file, "gro", output_gro_dir, output_gro_name
    )
    output_gro_file = export_gro_handler(gro_handler, output_gro_file, parser)
    paths = GromacsPaths(output_itp_file, output_gro_file, output_top_file)
    return paths


def process_solvent_itp(
    input_itp_file: str,
    new_residue_name: Optional[str] = None,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
):
    sections = parser.parse(input_itp_file)
    if new_residue_name:
        moleculetype_section = sections["data_moleculetype"]
        moleculetype_section = rename_data_column_content(
            moleculetype_section, "name", new_residue_name
        )
        atoms_section = sections["data_atoms"]
        atoms_section = rename_data_column_content(
            atoms_section, "res", new_residue_name
        )
        sections["data_moleculetype"] = moleculetype_section
        sections["data_atoms"] = atoms_section
    output_itp_path = prepare_output_file_path(
        input_itp_file, "itp", output_dir, output_name
    )
    output_itp_path = parser.export(sections, output_itp_path)
    return output_itp_path


def prepare_solvent_topol(
    input_top_file: str,
    residue_number: int,
    forcefield: str,
    new_include_file: str,
    new_residue_name: Optional[str] = None,  # New single include file for solvent
    output_name: Optional[str] = None,
    output_dir: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
    del_posre: bool = True,
    del_defaults: bool = True,
) -> str:
    residue_number = str(residue_number)
    sections = parser.parse(input_top_file)

    sections = delete_all_include_sections(sections)

    sections["include_solvent_itp"] = create_includes_section(new_include_file)
    sections.move_to_end("include_solvent_itp", last=False)

    sections["include_forcefield"] = create_includes_section(forcefield)
    sections.move_to_end("include_forcefield", last=False)

    if "data_molecules" not in sections:
        raise ValueError("No 'data_molecules' section found in topology file.")
    else:
        data_molecules_section = sections["data_molecules"]
        data_molecules_handler = parser.handler_registry.get_handler(
            data_molecules_section.construct_name
        )()
        data_molecules_handler.process(data_molecules_section)
        data_molecules_df = data_molecules_handler.content
        if len(data_molecules_df) != 1:
            raise ValueError("Multiple rows in 'data_molecules' section.")
        if new_residue_name:
            data_molecules_df["Compound"] = new_residue_name
        data_molecules_df["nmols"] = residue_number
        data_molecules_handler.content = data_molecules_df
        data_molecules_section = data_molecules_handler.export()
        sections["data_molecules"] = data_molecules_section

    # Optionally delete `data_posre` section
    if del_posre and "conditional_if" in sections:
        del sections["conditional_if"]

    if del_defaults and "data_defaults" in sections:
        del sections["data_defaults"]

    output_path = prepare_output_file_path(
        input_top_file, "top", output_dir, output_name
    )
    output_path = parser.export(sections, output_path)

    return output_path
