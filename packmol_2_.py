import logging
import subprocess
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)


# Helper functions
def validate_pdb_file(pdb_file: str) -> bool:
    """
    Validates a PDB file to ensure it has valid atom coordinates.
    """
    logger.info(f"Validating PDB file: {pdb_file}")
    try:
        with open(pdb_file, "r") as f:
            lines = f.readlines()
        # Check if at least one ATOM/HETATM line exists
        valid = any(line.startswith(("ATOM", "HETATM")) for line in lines)
        if not valid:
            logger.error("PDB file does not contain valid ATOM or HETATM entries.")
            return False
        logger.info("PDB file validation successful.")
        return True
    except FileNotFoundError:
        logger.error(f"PDB file not found: {pdb_file}")
        return False


def generate_packmol_input_for_solvent(
    solvent_file, output_file, num_molecules, box_size, tolerance=2.0
):
    """
    Generates a Packmol input script for solvent-only boxes.
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
    logger.debug(f"Packmol script:\n{packmol_script}")
    return script_path


def run_packmol(packmol_input: str, packmol_executable: str = "packmol") -> None:
    """
    Runs Packmol using the generated input script.
    """
    logger.info(f"Running Packmol with input file: {packmol_input}")
    if not Path(packmol_executable).exists():
        raise FileNotFoundError(f"Packmol binary not found: {packmol_executable}")

    try:
        with open(packmol_input, "r") as input_file:
            subprocess.run(
                [packmol_executable],
                stdin=input_file,
                text=True,
                capture_output=True,
                check=True,
            )
        logger.info("Packmol execution successful.")
    except subprocess.CalledProcessError as e:
        logger.error(f"Packmol failed:\n{e.stderr}")
        raise RuntimeError("Packmol execution failed.")


# Main test function
def test_packmol(
    solvent_file: str, output_file: str, box_size: list, num_molecules: int
):
    """
    Tests Packmol execution with a solvent file.
    """
    logger.info("Starting Packmol test...")

    # Step 1: Validate PDB file
    if not validate_pdb_file(solvent_file):
        raise ValueError("Invalid PDB file. Aborting Packmol test.")

    # Step 2: Generate Packmol input script
    packmol_input = generate_packmol_input_for_solvent(
        solvent_file=solvent_file,
        output_file=output_file,
        num_molecules=num_molecules,
        box_size=box_size,
    )

    # Step 3: Run Packmol
    try:
        run_packmol(packmol_input)
        logger.info(f"Packmol test completed successfully. Output: {output_file}")
    except RuntimeError:
        logger.error("Packmol test failed.")


# Test parameters
if __name__ == "__main__":
    solvent_file = "styrene.pdb"  # Replace with your solvent file
    output_file = "test_box.pdb"
    box_size = [30.0, 30.0, 30.0]  # Box size in Å
    num_molecules = 100

    try:
        test_packmol(solvent_file, output_file, box_size, num_molecules)
    except ValueError as e:
        logger.error(e)
