import os
import subprocess
import shutil
from typing import Dict
from processing.metadata_tracker import MetadataTracker
from data_models.solvent import Solvent


class SolvationBoxPreparer:
    """
    Prepares a simulation box using GROMACS solvation and ion addition, with metadata tracking.
    """

    ACCEPTABLE_DENSITY_TOLERANCE = 0.05  # ±5%

    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def prepare_box(self, solvent: Solvent, polymer_top_file: str, polymer_gro_file: str, box_size: float, output_dir: str, neutralize: bool = True) -> Dict[str, str]:
        """
        Prepare the simulation box by solvating with the specified solvent and optionally adding ions.

        Args:
            solvent (Solvent): Solvent object containing `.itp`, `.pdb`, and density/molecular weight.
            polymer_top_file (str): Path to the polymer topology file.
            polymer_gro_file (str): Path to the polymer `.gro` file (pre-simulation box).
            box_size (float): Size of the cubic box in nm.
            output_dir (str): Directory to save output files.
            neutralize (bool): Whether to neutralize the box with ions. Defaults to True.

        Returns:
            Dict[str, str]: Paths to generated files (`gro`, `top`, etc.).
        """
        os.makedirs(output_dir, exist_ok=True)

        polymer_name = os.path.splitext(os.path.basename(polymer_top_file))[0]

        # Backup the original topology file
        backup_top = os.path.join(output_dir, f"{polymer_name}_topology_backup.top")
        shutil.copy(polymer_top_file, backup_top)

        # Edit topology to include solvent itp and force field
        self._edit_topology(polymer_top_file, solvent)

        # Step 1: Create cubic box
        box_gro_file = os.path.join(output_dir, f"{polymer_name}_box.gro")
        self._create_cubic_box(polymer_gro_file, box_gro_file, box_size)

        # Step 2: Solvate the box with the solvent
        solvated_gro_file = os.path.join(output_dir, f"{polymer_name}_solvated.gro")
        self._solvate_box(box_gro_file, solvent.pdb_path, solvated_gro_file, polymer_top_file)

        # Step 3: Optionally neutralize the box
        if neutralize:
            neutralized_gro_file = os.path.join(output_dir, f"{polymer_name}_neutralized.gro")
            self._add_ions(solvated_gro_file, polymer_top_file, neutralized_gro_file)
        else:
            neutralized_gro_file = solvated_gro_file

        # Add metadata for the process
        self.metadata_tracker.add_step("Simulation Box Preparation", {
            "polymer_topology_file": polymer_top_file,
            "solvent_pdb_file": solvent.pdb_path,
            "box_gro_file": box_gro_file,
            "final_gro_file": neutralized_gro_file,
            "neutralized": neutralize
        })

        return {"final_gro_file": neutralized_gro_file, "topology_file": polymer_top_file}

    def _edit_topology(self, topology_file: str, solvent: Solvent):
        """
        Modify the topology file to include the solvent .itp and force field.

        Args:
            topology_file (str): Path to the topology file.
            solvent (Solvent): Solvent object containing necessary paths.
        """
        with open(topology_file, "a") as f:
            f.write(f'\n#include "{solvent.force_field}/forcefield.itp"\n')
            f.write(f'#include "{solvent.itp_path}"\n')

    def _create_cubic_box(self, input_gro: str, output_gro: str, box_size: float):
        """
        Create a cubic simulation box using GROMACS editconf.

        Args:
            input_gro (str): Path to the input .gro file.
            output_gro (str): Path to the output box .gro file.
            box_size (float): Desired box size in nm.
        """
        editconf_command = [
            "gmx", "editconf",
            "-f", input_gro,
            "-o", output_gro,
            "-c",  # Center the polymer
            "-d", str(box_size),  # Distance to the box edge
            "-bt", "cubic"  # Box type
        ]
        subprocess.run(editconf_command, check=True)

    def _solvate_box(self, box_gro: str, solvent_pdb: str, output_gro: str, topology_file: str):
        """
        Solvate the box using GROMACS solvate.

        Args:
            box_gro (str): Path to the cubic box .gro file.
            solvent_pdb (str): Path to the solvent .pdb file.
            output_gro (str): Path to the solvated .gro file.
            topology_file (str): Path to the topology file to update.
        """
        solvate_command = [
            "gmx", "solvate",
            "-cp", box_gro,
            "-cs", solvent_pdb,
            "-o", output_gro,
            "-p", topology_file
        ]
        subprocess.run(solvate_command, check=True)

    def _add_ions(self, input_gro: str, topology_file: str, output_gro: str):
        """
        Add ions to neutralize the simulation box using GROMACS genion.

        Args:
            input_gro (str): Path to the solvated .gro file.
            topology_file (str): Path to the topology file.
            output_gro (str): Path to the neutralized .gro file.
        """
        # Generate .tpr file for genion
        tpr_file = os.path.join(os.path.dirname(output_gro), "ions.tpr")
        grompp_command = [
            "gmx", "grompp",
            "-f", self._get_template("ions.mdp"),
            "-c", input_gro,
            "-p", topology_file,
            "-o", tpr_file
        ]
        subprocess.run(grompp_command, check=True)

        # Add Na+ and Cl− ions to neutralize
        genion_command = [
            "gmx", "genion",
            "-s", tpr_file,
            "-o", output_gro,
            "-p", topology_file,
            "-pname", "NA",
            "-nname", "CL",
            "-neutral"
        ]
        subprocess.run(genion_command, check=True)

    def _get_template(self, filename: str) -> str:
        """
        Retrieve the path to a template file (e.g., ions.mdp).

        Args:
            filename (str): Name of the template file.

        Returns:
            str: Full path to the template file.
        """
        template_dir = os.path.abspath("templates")
        template_path = os.path.join(template_dir, filename)
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file '{filename}' not found in {template_dir}.")
        return template_path
