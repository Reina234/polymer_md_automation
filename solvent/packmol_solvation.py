import subprocess
import os
from pdb_processing.managers.pdb_manager import PDBManager

class PackmolBoxGenerator:
    def __init__(self, pdb_manager: PDBManager, box_size: float, output_dir: str):
        """
        :param pdb_manager: Instance of PDBManager with a validated PDB.
        :param box_size: Size of the cubic simulation box (in Angstroms).
        :param output_dir: Directory to save the generated PDB file.
        """
        self.pdb_manager = pdb_manager
        self.box_size = box_size
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

    def calculate_molecule_count(self) -> int:
        """
        Calculate the number of solvent molecules required for the given density.
        """
        solvent = self.pdb_manager.solvent  # Get solvent metadata
        volume_cm3 = (self.box_size * 1e-8) ** 3  # Angstrom³ to cm³
        total_mass_g = solvent.density * volume_cm3  # Total mass in grams
        molecule_count = int((total_mass_g * 6.022e23) / solvent.molecular_weight)  # Avogadro's number
        return molecule_count

    def generate_box(self) -> str:
        """
        Generate the solvent box using Packmol.
        :return: Path to the generated Packmol PDB file.
        """
        molecule_count = self.calculate_molecule_count()
        input_pdb = self.pdb_manager.pdb_path # Get the validated PDB file
        output_pdb = os.path.join(self.output_dir, f"{self.pdb_manager.solvent.name}_box.pdb")

        packmol_input = f"""
        tolerance 2.0
        filetype pdb
        output {output_pdb}
        structure {input_pdb}
          number {molecule_count}
          inside box 0.0 0.0 0.0 {self.box_size} {self.box_size} {self.box_size}
        end structure
        """
        with open("packmol_input.inp", "w") as f:
            f.write(packmol_input)

        # Run Packmol
        result = subprocess.run(["packmol"], input="packmol_input.inp", text=True, capture_output=True)
        if result.returncode != 0:
            raise RuntimeError(f"Packmol failed:\n{result.stderr}")

        print(f"[+] Packmol box generated: {output_pdb}")
        return output_pdb
