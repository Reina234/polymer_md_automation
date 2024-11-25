# core/lammps_handler.py
import yaml
import os
import subprocess


class LAMMPSHandler:
    """
    Handles LAMMPS-specific tasks such as input script generation and simulation execution.
    """

    def __init__(self, config_file: str = "lammps_config.yaml"):
        """
        Initialize the handler with a configuration file.
        :param config_file: Path to the LAMMPS configuration file.
        """
        if not os.path.exists(config_file):
            raise FileNotFoundError(f"LAMMPS configuration file not found: {config_file}")
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)

    def write_input_script(self, data_file: str, output_file: str) -> None:
        """
        Write a LAMMPS input script using the provided configuration and data file.
        :param data_file: Path to the LAMMPS data file.
        :param output_file: Path to save the input script.
        """
        with open(output_file, "w") as f:
            f.write(f"read_data {data_file}\n")
            for key, value in self.config.items():
                if isinstance(value, list):  # Handle lists like "fix"
                    for line in value:
                        f.write(f"{key} {line}\n")
                else:
                    f.write(f"{key} {value}\n")
        print(f"LAMMPS input script written to {output_file}.")

    def run_simulation(self, input_file: str, lammps_exe: str = "lmp_mpi") -> None:
        """
        Run a LAMMPS simulation.
        :param input_file: Path to the LAMMPS input script.
        :param lammps_exe: Path to the LAMMPS executable.
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"LAMMPS input file not found: {input_file}")
        result = subprocess.run([lammps_exe, "-in", input_file], capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"LAMMPS simulation failed:\n{result.stderr}")
        print(f"LAMMPS simulation completed successfully.")
