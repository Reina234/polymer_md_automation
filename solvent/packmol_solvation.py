import subprocess
import os
from pdb_processing.managers.pdb_manager import PDBManager

#NOTE: change it so that input script is generated first, then running packmol checks if it exists, if not, it asks you to run generate


class PackmolBoxGenerator:
    def __init__(self, pdb_manager: PDBManager, box_size: float, output_dir: str):
        """
        :param pdb_manager: Instance of PDBManager with a validated PDB.
        :param box_size: Size of the cubic simulation box (in Angstroms).
        :param output_dir: Directory to save the generated PDB file.
        """
        self.pdb_manager = pdb_manager
        self.box_size = box_size
        self.output_dir = os.path.abspath(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)

    def _calculate_molecule_count(self) -> int:
        """
        Calculate the number of solvent molecules required for the given density.
        """
        solvent = self.pdb_manager.solvent  # Get solvent metadata
        volume_cm3 = (self.box_size * 1e-8) ** 3  # Convert Angstrom³ to cm³
        total_mass_g = solvent.density * volume_cm3  # Total mass in grams
        molecule_count = int((total_mass_g * 6.022e23) / solvent.molecular_weight)  # Avogadro's number
        
        return molecule_count

    def generate_box(self) -> str:
        """
        Generate the solvent box using Packmol.
        :return: Path to the generated Packmol PDB file.
        """
        molecule_count = self._calculate_molecule_count()
        input_pdb = self.pdb_manager.pdb_path  # Get the validated PDB file
        output_pdb = os.path.join(self.output_dir, f"{self.pdb_manager.solvent.name}_box.pdb")

        # Generate the Packmol input file
        packmol_input = f"""
        tolerance 2.0
        filetype pdb
        output {output_pdb}
        structure {input_pdb}
          number {molecule_count}
          inside box 0.0 0.0 0.0 {self.box_size} {self.box_size} {self.box_size}
        end structure
        """
        input_file_path = os.path.join(self.output_dir, "packmol_input.inp")
        with open(input_file_path, "w", newline="\n") as f:
            f.write(packmol_input)

        command = f"packmol < {input_file_path}"
        exit_code = os.system(command)
        if exit_code == 0:
            print(f"[INFO] Packmol executed successfully. File saved at {output_pdb}.")
        else:
            print(f"[ERROR] Packmol execution failed with exit code {exit_code}.")
