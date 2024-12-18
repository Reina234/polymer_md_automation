from typing import List
import subprocess
import logging
from A_modules.shared.utils.calculation_utils import calculate_num_particles
from A_config.constants import LengthUnits2, MassUnits2

logger = logging.getLogger(__name__)


def create_density_accurate_solvent_box(
    solvent_file: str,
    output_file: str,
    box_size: List[float],
    target_density: float,
    molecular_weight: float,
    tolerance: float = 2.0,  # Å
    packmol_executable: str = "packmol",
    editconf_executable: str = "gmx",
) -> str:
    """
    Creates a density-accurate solvent-only box using Packmol.

    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path for the generated solvent box file (.gro).
    :param box_size: Dimensions of the box [x, y, z] in Å.
    :param target_density: Desired system density in g/cm³.
    :param molecular_weight: Molecular weight of the solvent in g/mol.
    :param tolerance: Minimum distance tolerance in Å (default: 2.0).
    :param packmol_executable: Path to the Packmol binary (default: "packmol").
    :param editconf_executable: Path to the GROMACS binary (default: "gmx").
    :return: Path to the generated .gro file.
    """
    # Step 1: Calculate the number of solvent molecules required
    target_num_molecules = calculate_num_particles(
        box_dimensions=box_size,
        molecular_weight=molecular_weight,
        density_SI=target_density,
        box_units=LengthUnits2.ANGSTROM,
        mass_units=MassUnits2.GRAM,
    )

    logger.info(f"Target number of solvent molecules: {target_num_molecules}")

    # Step 2: Generate Packmol input script
    packmol_input = generate_packmol_input_for_solvent(
        solvent_file=solvent_file,
        output_file=output_file.replace(".gro", ".pdb"),
        num_molecules=target_num_molecules,
        box_size=box_size,
        tolerance=tolerance,
    )

    # Step 3: Run Packmol
    run_packmol(packmol_input, packmol_executable)

    # Step 4: Convert .pdb to .gro
    gro_output_file = output_file
    pdb_output_file = output_file.replace(".gro", ".pdb")
    convert_pdb_to_gro(
        input_pdb=pdb_output_file,
        output_gro=gro_output_file,
        editconf_executable=editconf_executable,
    )

    logger.info(f"Solvent box created and saved to: {gro_output_file}")
    return gro_output_file


def generate_packmol_input_for_solvent(
    solvent_file: str,
    output_file: str,
    num_molecules: int,
    box_size: List[float],
    tolerance: float = 2.0,
) -> str:
    """
    Generates a Packmol input script for a solvent-only box.

    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path for the output solvated box file (.pdb).
    :param num_molecules: Number of solvent molecules to place in the box.
    :param box_size: Dimensions of the box [x, y, z] in Å.
    :param tolerance: Minimum distance tolerance in Å (default: 2.0).
    :return: Path to the generated Packmol input script.
    """
    packmol_script = f"""tolerance {tolerance}
filetype pdb

output {output_file}

structure {solvent_file}
  number {num_molecules}
  inside box 0.0 0.0 0.0 {box_size[0]} {box_size[1]} {box_size[2]}
end structure
"""

    script_path = output_file.replace(".pdb", ".inp")
    with open(script_path, "w") as f:
        f.write(packmol_script)

    return script_path


def generate_packmol_input_for_solvent(
    solvent_file: str,
    output_file: str,
    num_molecules: int,
    box_size: List[float],
    tolerance: float = 2.0,
) -> str:
    """
    Generates a Packmol input script for a solvent-only box.

    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path for the output solvated box file (.pdb).
    :param num_molecules: Number of solvent molecules to place in the box.
    :param box_size: Dimensions of the box [x, y, z] in Å.
    :param tolerance: Minimum distance tolerance in Å (default: 2.0).
    :return: Path to the generated Packmol input script.
    """
    packmol_script = f"""tolerance {tolerance}
filetype pdb

output {output_file}

structure {solvent_file}
  number {num_molecules}
  inside box 0.0 0.0 0.0 {box_size[0]} {box_size[1]} {box_size[2]}
end structure
"""

    script_path = output_file.replace(".pdb", ".inp")
    with open(script_path, "w") as f:
        f.write(packmol_script)

    return script_path


def generate_packmol_input_for_solvent(
    solvent_file: str,
    output_file: str,
    num_molecules: int,
    box_size: List[float],
    tolerance: float = 2.0,
) -> str:
    """
    Generates a Packmol input script for a solvent-only box.

    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path for the output solvated box file (.pdb).
    :param num_molecules: Number of solvent molecules to place in the box.
    :param box_size: Dimensions of the box [x, y, z] in Å.
    :param tolerance: Minimum distance tolerance in Å (default: 2.0).
    :return: Path to the generated Packmol input script.
    """
    packmol_script = f"""tolerance {tolerance}
filetype pdb

output {output_file}

structure {solvent_file}
  number {num_molecules}
  inside box 0.0 0.0 0.0 {box_size[0]} {box_size[1]} {box_size[2]}
end structure
"""

    script_path = output_file.replace(".pdb", ".inp")
    with open(script_path, "w") as f:
        f.write(packmol_script)

    return script_path


def run_packmol(packmol_input: str, packmol_executable: str = "packmol") -> None:
    """
    Runs Packmol to generate a solvated box.
    """
    try:
        result = subprocess.run(
            [packmol_executable, "<", packmol_input],
            check=True,
            text=True,
            capture_output=True,
            shell=True,
        )
        logger.info(result.stdout)
    except subprocess.CalledProcessError as e:
        logger.error(f"Packmol failed: {e.stderr}")
        raise RuntimeError("Packmol execution failed.")


def convert_pdb_to_gro(
    input_pdb: str, output_gro: str, editconf_executable: str = "gmx"
) -> None:
    """
    Converts a .pdb file to .gro format using GROMACS editconf.

    :param input_pdb: Path to the input .pdb file.
    :param output_gro: Path to the output .gro file.
    :param editconf_executable: Path to the GROMACS binary (default: "gmx").
    """
    try:
        subprocess.run(
            [
                editconf_executable,
                "editconf",
                "-f",
                input_pdb,
                "-o",
                output_gro,
            ],
            check=True,
        )
        logger.info(f"Converted {input_pdb} to {output_gro}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Conversion to .gro failed: {e.stderr}")
        raise RuntimeError("GROMACS editconf execution failed.")
