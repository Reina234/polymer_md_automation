from modules.shared.utils.file_utils import (
    file_type_check_wrapper,
    prepare_output_file_path,
)
from collections import OrderedDict
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler
from typing import Optional, Tuple
import pandas as pd
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
import logging
import os
import shutil
from data_models.output_types import GromacsPaths
from modules.atomistic.utils.file_utils import (
    get_gro_handler,
    get_residue_number,
    rename_residue_name_from_handler,
    rename_data_column_content,
    calculate_molecule_counts,
    replace_dataframe_contents,
    export_gro_handler,
    rename_residue_name_from_gro,
    replace_value_in_dataframe,
    create_includes_section,
    delete_all_include_sections,
    validate_and_extract_residue_name,
    add_full_rows_to_handler_deduplicate,
    add_to_specific_handler_columns,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def prepare_solute_files(
    solute_itp_file: str,
    solvent_itp_file: str,
    solvent_box_gro_file: str,
    input_top_file: str,
    output_dir: str,
    forcefield: str = "amber99sb-ildn.ff/forcefield.itp",
    ion_itp_file: str = "amber99sb-ildn.ff/ions.itp",
    solute_molecule_name: str = "POLY",
    output_solvent_gro_name: Optional[str] = None,
    output_solute_itp_name: Optional[str] = None,
    output_solvent_itp_name: Optional[str] = None,
    output_top_name: str = "topol",
    parser: GromacsParser = GromacsParser(),
) -> GromacsPaths:
    print(solvent_box_gro_file)
    gro_handler = get_gro_handler(solvent_box_gro_file)
    # NOTE: need to fix to allow for separate things
    # residue_number = get_residue_number(gro_handler)
    # residue_name = validate_and_extract_residue_name(gro_handler)
    molecule_content = calculate_molecule_counts(
        gro_handler=gro_handler,
        residue_name_col="Residue Name",
        residue_number_col="Residue Number",
    )
    molecule_content = replace_value_in_dataframe(
        molecule_content,
        target_value="UNL",
        replacement_value=solute_molecule_name,
        move_to_top=True,
    )

    solute_itp_file, solvent_itp_file = process_solute_and_solvent_itps(
        solute_itp_file=solute_itp_file,
        solvent_itp_file=solvent_itp_file,
        output_dir=output_dir,
        solute_output_name=output_solute_itp_name,
        solvent_output_name=output_solvent_itp_name,
        parser=parser,
    )

    solute_itp_path = os.path.abspath(solute_itp_file)
    solvent_itp_path = os.path.abspath(solvent_itp_file)

    output_top_path = prepare_solute_topol(
        input_top_file=input_top_file,
        new_molecule_dataframe=molecule_content,
        forcefield=forcefield,
        ions_itp_file=ion_itp_file,
        solute_itp_file=solute_itp_path,
        solvent_itp_file=solvent_itp_path,
        output_name=output_top_name,
        output_dir=output_dir,
        parser=parser,
        del_posre=True,
        del_defaults=True,
    )

    output_gro_path = prepare_output_file_path(
        solvent_box_gro_file, output_extension="gro", output_dir=output_dir
    )
    if output_gro_path != solvent_box_gro_file:
        shutil.copy(solvent_box_gro_file, output_gro_path)
    paths = GromacsPaths(
        itp_path=solute_itp_file, gro_path=output_gro_path, top_path=output_top_path
    )
    return paths


def process_solvent_itp(
    solvent_itp_file: str,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
    data_handler: DataHandler = DataHandler,
) -> Tuple[str, pd.DataFrame]:
    data_handler = data_handler()

    solvent_sections = parser.parse(solvent_itp_file)
    atoms_sections = solvent_sections["data_atomtypes"]
    data_handler.process(atoms_sections)
    atom_content = data_handler.content
    solvent_sections.pop("data_atomtypes")
    output_itp_path = prepare_output_file_path(
        solvent_itp_file, "itp", output_dir, output_name
    )
    output_itp_path = parser.export(solvent_sections, output_itp_path)
    return output_itp_path, atom_content


def process_solute_itp(
    solute_itp_file: str,
    solvent_atomtype_data: pd.DataFrame,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
    data_handler: DataHandler = DataHandler,
) -> str:
    data_handler = data_handler()
    solute_sections = parser.parse(solute_itp_file)
    atoms_sections = solute_sections["data_atomtypes"]

    data_handler.process(atoms_sections)
    data_handler = add_full_rows_to_handler_deduplicate(
        data_handler, solvent_atomtype_data, add_to_top=False, deduplicate_column="name"
    )
    atoms_sections = data_handler.export()
    solute_sections["data_atomtypes"] = atoms_sections

    output_itp_path = prepare_output_file_path(
        solute_itp_file, "itp", output_dir, output_name
    )
    output_itp_path = parser.export(solute_sections, output_itp_path)
    return output_itp_path


def process_solute_and_solvent_itps(
    solute_itp_file: str,
    solvent_itp_file: str,
    output_dir: str,
    solute_output_name: str,
    solvent_output_name: str,
    parser: GromacsParser = GromacsParser(),
    data_handler: DataHandler = DataHandler,
) -> Tuple[str, str]:
    output_solvent_itp, solvent_atom_data = process_solvent_itp(
        solvent_itp_file, output_dir, solvent_output_name, parser=parser
    )
    output_solute_itp = process_solute_itp(
        solute_itp_file,
        solvent_atom_data,
        output_dir=output_dir,
        output_name=solute_output_name,
        parser=parser,
    )

    return output_solute_itp, output_solvent_itp


def prepare_solute_topol(
    input_top_file: str,
    new_molecule_dataframe: pd.DataFrame,
    forcefield: str,
    ions_itp_file: str,
    solute_itp_file: str,
    solvent_itp_file: str,
    output_name: Optional[str] = None,
    output_dir: Optional[str] = None,
    parser: GromacsParser = GromacsParser(),
    del_posre: bool = True,
    del_defaults: bool = True,
):

    sections = parser.parse(input_top_file)

    sections = delete_all_include_sections(sections)

    sections["include_ions_itp"] = create_includes_section(ions_itp_file)
    sections.move_to_end("include_ions_itp", last=False)

    sections["include_solvent_itp"] = create_includes_section(solvent_itp_file)
    sections.move_to_end("include_solvent_itp", last=False)

    sections["include_solute_itp"] = create_includes_section(solute_itp_file)
    sections.move_to_end("include_solute_itp", last=False)

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
        new_content = replace_dataframe_contents(
            original_df=data_molecules_handler.content,
            new_df=new_molecule_dataframe,
            pad_missing=True,
        )
        data_molecules_handler.content = new_content
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
