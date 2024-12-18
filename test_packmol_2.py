import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_pdb_file(pdb_file: str) -> None:
    """
    Validate the PDB file to ensure it has valid ATOM or HETATM entries.
    :param pdb_file: Path to the PDB file.
    :raises ValueError: If the PDB file is invalid.
    """
    logger.info(f"Validating PDB file: {pdb_file}")
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    valid = any(line.startswith(("ATOM", "HETATM")) for line in lines)
    if not valid:
        raise ValueError(
            f"PDB file {pdb_file} does not contain valid ATOM or HETATM entries."
        )
    logger.info(f"PDB file {pdb_file} is valid.")


def calculate_molecule_count(density: float, molar_mass: float, box_size: float) -> int:
    """
    Calculate the number of molecules required to achieve the target density.
    :param density: Target density in g/cm³.
    :param molar_mass: Molar mass of the solvent in g/mol.
    :param box_size: Length of the cubic simulation box in Angstroms.
    :return: Number of molecules.
    """
    volume_cm3 = (box_size * 1e-8) ** 3  # Convert Å³ to cm³
    total_mass_g = density * volume_cm3  # Total mass in grams
    num_molecules = int((total_mass_g * 6.022e23) / molar_mass)  # Avogadro's number
    logger.info(
        f"Calculated {num_molecules} molecules for box size {box_size} Å³ and density {density} g/cm³."
    )
    return num_molecules


def generate_packmol_input(
    solvent_file: str, output_file: str, num_molecules: int, box_size: float
) -> str:
    """
    Generate the Packmol input script.
    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path to the output PDB file.
    :param num_molecules: Number of molecules to place in the box.
    :param box_size: Size of the cubic box in Angstroms.
    :return: Path to the generated Packmol input script.
    """
    packmol_input = f"""
tolerance 2.0
filetype pdb

output {output_file}

structure {solvent_file}
  number {num_molecules}
  inside box 0.0 0.0 0.0 {box_size} {box_size} {box_size}
end structure
"""
    input_file = output_file.replace(".pdb", ".inp")
    with open(input_file, "w", newline="\n") as f:
        f.write(packmol_input)
    logger.info(f"Generated Packmol input script:\n{packmol_input}")
    return input_file


def run_packmol(input_script: str) -> None:
    """
    Run Packmol using the generated input script.
    :param input_script: Path to the Packmol input script.
    """
    logger.info(f"Running Packmol with input file: {input_script}")
    if not os.path.isfile(input_script):
        raise FileNotFoundError(f"Packmol input script not found: {input_script}")

    command = f"packmol < {input_script}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Packmol failed with return code {result.returncode}.")
        logger.error(f"Packmol stderr:\n{result.stderr}")
        raise RuntimeError("Packmol execution failed.")
    logger.info(f"Packmol executed successfully. Output saved.")


import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_pdb_file(pdb_file: str) -> None:
    """
    Validate the PDB file to ensure it has valid ATOM or HETATM entries.
    :param pdb_file: Path to the PDB file.
    :raises ValueError: If the PDB file is invalid.
    """
    logger.info(f"Validating PDB file: {pdb_file}")
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    with open(pdb_file, "r") as f:
        lines = f.readlines()

    valid = any(line.startswith(("ATOM", "HETATM")) for line in lines)
    if not valid:
        raise ValueError(
            f"PDB file {pdb_file} does not contain valid ATOM or HETATM entries."
        )
    logger.info(f"PDB file {pdb_file} is valid.")


def calculate_molecule_count(density: float, molar_mass: float, box_size: float) -> int:
    """
    Calculate the number of molecules required to achieve the target density.
    :param density: Target density in g/cm³.
    :param molar_mass: Molar mass of the solvent in g/mol.
    :param box_size: Length of the cubic simulation box in Angstroms.
    :return: Number of molecules.
    """
    volume_cm3 = (box_size * 1e-8) ** 3  # Convert Å³ to cm³
    total_mass_g = density * volume_cm3  # Total mass in grams
    num_molecules = int((total_mass_g * 6.022e23) / molar_mass)  # Avogadro's number
    logger.info(
        f"Calculated {num_molecules} molecules for box size {box_size} Å³ and density {density} g/cm³."
    )
    return num_molecules


def generate_packmol_input(
    solvent_file: str, output_file: str, num_molecules: int, box_size: float
) -> str:
    """
    Generate the Packmol input script.
    :param solvent_file: Path to the solvent PDB file.
    :param output_file: Path to the output PDB file.
    :param num_molecules: Number of molecules to place in the box.
    :param box_size: Size of the cubic box in Angstroms.
    :return: Path to the generated Packmol input script.
    """
    packmol_input = f"""
tolerance 2.0
filetype pdb

output {output_file}

structure {solvent_file}
  number {num_molecules}
  inside box 0.0 0.0 0.0 {box_size} {box_size} {box_size}
end structure
"""
    input_file = output_file.replace(".pdb", ".inp")
    with open(input_file, "w", newline="\n") as f:
        f.write(packmol_input)
    logger.info(f"Generated Packmol input script:\n{packmol_input}")
    return input_file


def run_packmol(input_script: str) -> None:
    """
    Run Packmol using the generated input script.
    :param input_script: Path to the Packmol input script.
    """
    logger.info(f"Running Packmol with input file: {input_script}")
    if not os.path.isfile(input_script):
        raise FileNotFoundError(f"Packmol input script not found: {input_script}")

    command = f"packmol < {input_script}"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Packmol failed with return code {result.returncode}.")
        logger.error(f"Packmol stderr:\n{result.stderr}")
        raise RuntimeError("Packmol execution failed.")
    logger.info(f"Packmol executed successfully. Output saved.")


### Test Script
def convert_pdb_to_gro(pdb_file: str, gro_file: str) -> None:
    """
    Convert a PDB file to a GRO file using GROMACS's editconf tool.
    :param pdb_file: Path to the input PDB file.
    :param gro_file: Path to the output GRO file.
    """
    logger.info(f"Converting {pdb_file} to {gro_file} using GROMACS editconf.")
    if not os.path.isfile(pdb_file):
        raise FileNotFoundError(f"PDB file not found: {pdb_file}")

    command = ["gmx", "editconf", "-f", pdb_file, "-o", gro_file]
    result = subprocess.run(command, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"GROMACS editconf failed with return code {result.returncode}.")
        logger.error(f"editconf stderr:\n{result.stderr}")
        raise RuntimeError("GROMACS editconf execution failed.")
    logger.info(f"Conversion successful. GRO file saved at {gro_file}.")


def test_generate_solvent_box():
    """
    Test the solvent box generation process with Packmol.
    """
    # Test parameters
    solvent_file = "styrene.pdb"  # Replace with your solvent PDB
    output_file = "test_box.pdb"
    box_size = 30.0  # Box size in Å
    density = 0.9  # Density in g/cm³ (example for styrene)
    molar_mass = 104.15  # Molar mass in g/mol (example for styrene)

    try:
        # Validate the PDB file
        validate_pdb_file(solvent_file)

        # Calculate the number of molecules needed
        num_molecules = calculate_molecule_count(density, molar_mass, box_size)

        # Generate the Packmol input script
        packmol_input = generate_packmol_input(
            solvent_file, output_file, num_molecules, box_size
        )

        # Run Packmol to generate the box
        run_packmol(packmol_input)

        logger.info(
            f"Test completed successfully. Generated box saved at {output_file}."
        )

    except (ValueError, FileNotFoundError, RuntimeError) as e:
        logger.error(f"Test failed: {e}")


if __name__ == "__main__":
    test_generate_solvent_box()
