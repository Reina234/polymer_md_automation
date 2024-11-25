import os
import yaml
from core.packmol_box_generator import PackmolBoxGenerator
from core.lammps_equilibrator import LAMMPSEquilibrator


class EquilibrationWorkflow:
    """
    Orchestrates solvent preparation, box generation, and LAMMPS equilibration.
    """

    def __init__(self, solvent, config_file: str):
        """
        :param solvent: Solvent object containing metadata and PDB file.
        :param config_file: Path to the YAML configuration file.
        """
        self.solvent = solvent
        with open(config_file, "r") as f:
            self.config = yaml.safe_load(f)
        self.output_dir = self.config["output_directory"]
        self.lammps_params = self.config["lammps_parameters"]
        os.makedirs(self.output_dir, exist_ok=True)

    def prepare(self, box_size: float, topology_file: str) -> str:
        """
        Prepare the solvent box and convert it to LAMMPS data format.
        :param box_size: Size of the cubic simulation box (in Angstroms).
        :param topology_file: Path to the topology/force field file.
        :return: Path to the generated LAMMPS data file.
        """
        print("Generating solvent box...")
        box_generator = PackmolBoxGenerator(self.solvent, box_size, self.output_dir)
        packmol_pdb = box_generator.generate_box()

        print("Converting to LAMMPS data format...")
        data_file = os.path.join(self.output_dir, f"{self.solvent.name}_data.lmp")
        self.solvent.generate_lammps_data(topology_file, data_file)

        return data_file

    def equilibrate(self, data_file: str) -> None:
        """
        Run the equilibration in LAMMPS.
        :param data_file: Path to the LAMMPS data file.
        """
        equilibrator = LAMMPSEquilibrator(self.lammps_params, self.output_dir)
        input_script = equilibrator.write_input_script(data_file)
        equilibrator.run_equilibration(input_script)
