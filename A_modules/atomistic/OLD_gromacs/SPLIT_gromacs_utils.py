from config.paths import TOPOL_NAME, BASE_OUTPUT_DIR, GROMACS_OUTPUT_SUBDIR
import os
from pathlib import Path
from typing import Optional
from preprocessing.parsers.top_parser import TOPParser
from preprocessing.parsers.itp_parser import ITPParser
from preprocessing.utils import copy_files
from preprocessing.solvent_atomtypes_manager import AtomtypesManager
import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def prepare_topol_file(input_top_path: str, run_name: str):
    output_dir = os.path.join(run_name, GROMACS_OUTPUT_SUBDIR)
    os.makedirs(output_dir, exist_ok=True)
    topol_file = move_and_rename_topol_file(input_top_path, output_dir, TOPOL_NAME)
    return topol_file


def move_and_rename_topol_file(
    input_top_path: str, output_dir: str, output_top_name: str
) -> str:
    """
    Copies and renames the topology file to the specified directory.

    Args:
        input_top_path (str): Path to the source topology file.
        output_dir (str): Path to the output directory.
        output_top_name (str): New name for the topology file.

    Returns:
        str: Full path to the renamed topology file.
    """
    # Ensure the source file exists
    if not os.path.exists(input_top_path):
        raise FileNotFoundError(
            f"The input topology file '{input_top_path}' does not exist."
        )

    top_name = os.path.basename(input_top_path)
    top_directory = os.path.dirname(input_top_path)

    # Copy the file to the output directory
    copied_files = copy_files([top_name], top_directory, output_dir)

    # Ensure at least one file is copied
    if not copied_files:
        raise ValueError(f"File '{top_name}' could not be copied to '{output_dir}'.")

    copied_file_path = Path(
        copied_files[0]
    )  # Extract the first file path from the list

    # Skip renaming if the source and destination paths are the same
    dest_file = os.path.join(output_dir, output_top_name)
    if os.path.abspath(copied_file_path) == os.path.abspath(dest_file):
        logger.info(
            f"Skipping rename: '{copied_file_path}' and '{dest_file}' are the same file."
        )
        return str(dest_file)

    # Rename the copied file
    renamed_file = copied_file_path.rename(dest_file)
    logger.info(f"Renamed '{copied_file_path}' to '{renamed_file}'.")
    return str(renamed_file)


def reformat_topol_file(
    input_top_path: str,
    solute_itp_path: str,
    forcefield_path: str,
    output_dir: Optional[str] = None,
    solvent_itp_path: Optional[str] = None,
    posres: bool = False,
) -> str:
    """
    Reformat the topology file to ensure correct include order and handle position restraints.

    Args:
        input_top_path (str): Path to the input topology file.
        solute_itp_path (str): Path to the solute .itp file.
        forcefield_path (str): Path to the forcefield .itp file.
        output_dir (Optional[str]): Directory to save the reformatted topology file.
        solvent_itp_path (Optional[str]): Path to the solvent .itp file (if used).
        posres (bool): Whether to include position restraints.

    Returns:
        str: Path to the reformatted topology file.
    """
    # Resolve absolute paths
    solute_itp_path = os.path.abspath(solute_itp_path)
    if solvent_itp_path:
        solvent_itp_path = os.path.abspath(solvent_itp_path)

    parser = TOPParser()

    # Ensure the input topology file is read
    content = parser.read_file(input_top_path)

    # Reformat includes and remove the defaults section
    content = parser.ensure_include_order(
        content,
        forcefield_include=forcefield_path,
        monomer_include=solute_itp_path,
        solvent_include=solvent_itp_path,
    )
    content = parser.handle_posres(content, posres)
    content = parser.remove_section(content, "defaults")

    # Save the updated topology file
    output_top_path = (
        os.path.join(output_dir, "topol.top") if output_dir else input_top_path
    )
    parser.save(output_top_path, content)
    return output_top_path


def add_atomtypes_to_topology(
    solvent_name: str, target_topology: str, include_header: bool = False
):
    """
    Add solvent atomtypes directly into the [ atomtypes ] section of the target topology file,
    while preserving any existing atomtypes.

    Args:
        solvent_name (str): Name of the solvent whose atomtypes to retrieve.
        target_topology (str): Path to the target topology file.
        include_header (bool): Whether to include the `[ atomtypes ]` header and comments.
    """
    manager = AtomtypesManager()
    new_atomtypes = manager.retrieve_atomtypes(solvent_name)

    # Optionally remove the `[ atomtypes ]` header and comments
    if not include_header:
        new_atomtypes = [
            line
            for line in new_atomtypes
            if not line.strip().startswith("[") and not line.strip().startswith(";")
        ]

    with open(target_topology, "r") as file:
        content = file.readlines()

    in_atomtypes_section = False
    updated_content = []
    existing_atomtypes = []
    atomtypes_added = False

    for line in content:
        stripped_line = line.strip()

        # Detect the start of the [ atomtypes ] section
        if stripped_line.startswith("[ atomtypes ]"):
            in_atomtypes_section = True
            updated_content.append(line)  # Keep the section header
            continue

        # Collect existing atomtypes within the section
        if in_atomtypes_section:
            if stripped_line.startswith("[") and not stripped_line.startswith(
                "[ atomtypes ]"
            ):
                in_atomtypes_section = False  # End of section
                updated_content.extend(existing_atomtypes)  # Add merged atomtypes
                updated_content.extend(new_atomtypes)  # Add new atomtypes
                atomtypes_added = True
                updated_content.append("\n")  # Add a blank line for separation
                updated_content.append(line)  # Keep the next section header
            else:
                existing_atomtypes.append(line)
                continue

        # Add lines to the updated content
        if not in_atomtypes_section:
            updated_content.append(line)

    # If no [ atomtypes ] section was found, create one at the beginning
    if not atomtypes_added:
        updated_content.insert(0, "\n")
        updated_content.insert(0, "\n".join(new_atomtypes) + "\n")
        updated_content.insert(0, "\n".join(existing_atomtypes) + "\n")
        updated_content.insert(0, "[ atomtypes ]\n")
        logger.info("[+] Created new [ atomtypes ] section in the topology file.")

    # Save the updated topology file
    with open(target_topology, "w") as file:
        file.writelines(updated_content)

    logger.info(
        f"[+] Added {solvent_name} atomtypes to the [ atomtypes ] section of {target_topology}."
    )
    return target_topology


def reformat_topol_file_for_solvent(
    input_top_path: str,
    forcefield_path: str,
    solvent_itp_path: str,
    solvent_name: str,
    num_molecules: int,
    output_dir: Optional[str] = None,
) -> str:
    """
    Reformat the topology file for solvent equilibrium by updating the number of solvent molecules.

    Args:
        input_top_path (str): Path to the input topology file.
        forcefield_path (str): Path to the forcefield .itp file.
        solvent_itp_path (str): Path to the solvent .itp file.
        solvent_name (str): Name of the solvent as specified in the .itp file.
        num_molecules (int): Number of solvent molecules to add.
        output_dir (Optional[str]): Directory to save the reformatted topology file.

    Returns:
        str: Path to the reformatted topology file.
    """
    parser = TOPParser()

    # Read the topology file
    content = parser.read_file(input_top_path)

    # Ensure includes are in the correct order
    content = parser.ensure_include_order(
        content,
        forcefield_include=forcefield_path,
        monomer_include=None,  # No solute for solvent equilibrium
        solvent_include=solvent_itp_path,
    )

    # Update the [ molecules ] section
    in_molecules_section = False
    updated_content = []
    molecule_updated = False

    for line in content:
        stripped_line = line.strip()

        if stripped_line.startswith("[ molecules ]"):
            in_molecules_section = True
            updated_content.append(line)
            continue

        if in_molecules_section:
            if stripped_line and not stripped_line.startswith(";"):
                molecule_data = stripped_line.split()
                if molecule_data[0] == solvent_name:
                    updated_content.append(f"{solvent_name} {num_molecules}\n")
                    molecule_updated = True
                else:
                    updated_content.append(line)
            else:
                in_molecules_section = False

        if not in_molecules_section:
            updated_content.append(line)

    # If solvent not found, add it to the [ molecules ] section
    if not molecule_updated:
        updated_content.append("\n[ molecules ]\n")
        updated_content.append(f"{solvent_name} {num_molecules}\n")

    # Save the updated topology file
    output_top_path = (
        os.path.join(output_dir, "topol.top") if output_dir else input_top_path
    )
    parser.save(output_top_path, updated_content)
    return output_top_path


from config.paths import MDP_DIRS, MDP_NAMING_SCHEME, MDP_TEMPLATE_PATHS, TemplatedMdps
from typing import Optional
import os
import logging


def create_mdps(
    mdp_type: TemplatedMdps,
    simulation_temp_k: float,
    template_path: Optional[str] = None,
    output_dir: Optional[str] = None,
):
    if not template_path:
        template_path = MDP_TEMPLATE_PATHS[mdp_type.value]

    if not output_dir:
        output_dir = MDP_DIRS[mdp_type.value]

    if not os.path.exists(template_path):
        raise FileNotFoundError(f"Template path {template_path} does not exist")

    # Construct the output file path
    output_filename = MDP_NAMING_SCHEME.format(
        run_type=mdp_type.value, temp=simulation_temp_k
    )
    output_path = os.path.join(output_dir, output_filename)

    # Check if the file already exists
    if os.path.exists(output_path):
        logger.info(f"MDP file already exists at {output_path}, skipping creation.")
        return output_path

    # Read and modify the template content
    with open(template_path, "r") as file:
        content = file.read()

    content = content.replace("{temp}", str(simulation_temp_k))
    os.makedirs(output_dir, exist_ok=True)

    # Write the modified content to the output file
    with open(output_path, "w") as file:
        file.write(content)

    logger.info(f"Created MDP file at {output_path}")
    return output_path


def retrieve_mdps(
    mdp_type: TemplatedMdps,
    simulation_temp_k: float,
    searched_dir: Optional[str] = None,
):
    if not searched_dir:
        searched_dir = MDP_DIRS[mdp_type.value]

    output_filename = MDP_NAMING_SCHEME.format(
        run_type=mdp_type.value, temp=simulation_temp_k
    )
    output_path = os.path.join(searched_dir, output_filename)

    if not os.path.exists(output_path):
        logger.warning(
            f"MDP file {output_path} does not exist, creating a new one instead."
        )
        return create_mdps(mdp_type, simulation_temp_k)

    logger.info(f"Retrieved MDP file from {output_path}")
    return output_path
