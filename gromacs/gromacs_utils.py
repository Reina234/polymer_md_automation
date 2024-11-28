from config.paths import TOPOL_NAME, BASE_OUTPUT_DIR, GROMACS_OUTPUT_SUBDIR
import os
from pathlib import Path
from typing import Optional
from preprocessing.parsers.top_parser import TOPParser
from preprocessing.utils import copy_files


def prepare_topol_file(
    input_top_path: str, run_name: str, output_base_dir: str = BASE_OUTPUT_DIR
):
    output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
    os.makedirs(output_dir, exist_ok=True)
    topol_file = move_and_rename_topol_file(input_top_path, output_dir, TOPOL_NAME)
    return topol_file


def move_and_rename_topol_file(
    input_top_path: str, output_dir: str, output_top_name: str
) -> str:
    top_name = os.path.basename(input_top_path)
    top_directory = os.path.dirname(input_top_path)
    topol_file = Path(copy_files([top_name], top_directory, output_dir))
    topol_file = topol_file.rename(os.path.join(output_dir, output_top_name))
    return str(topol_file)


def reformat_topol_file(
    input_top_path: str,
    solute_itp_path: str,
    solvent_itp_path: str,
    forcefield_path: str,
    output_dir: Optional[str] = None,
    posres: bool = False,
) -> str:
    parser = TOPParser()
    forcefield_include = forcefield_path
    solvent_include = solvent_itp_path
    solute_include = solute_itp_path
    content = parser.read_file(input_top_path)
    content = parser.ensure_include_order(
        content, forcefield_include, solute_include, solvent_include
    )
    parser.handle_posres(content, posres)
    content = parser.remove_section(content, "defaults")

    output_top_path = (
        os.path.join(output_dir, TOPOL_NAME) if output_dir else input_top_path
    )

    parser.save(output_top_path, content)
    return output_top_path
