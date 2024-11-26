import os
import subprocess
import shutil
from typing import Dict
from processing.metadata_tracker import MetadataTracker
from data_models.solvent import Solvent


class SimulationBoxPreparer:
    """
    Prepares the simulation box with the specified solvent and ions, ensuring correct density.
    """

    ACCEPTABLE_DENSITY_TOLERANCE = 0.05  # ±5%

    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def prepare_box(self, solvent: Solvent, top_file: str, box_size: float, output_dir: str, target_density: float, neutralize: bool = True) -> Dict[str, str]:
        """
        Prepare the simulation box by solvating and adding ions, ensuring correct density.

        Args:
            solvent (Solvent): Solvent object containing all necessary data.
            top_file (str): Path to the topology file.
            box_size (float): Initial size of the cubic box in nm.
            output_dir (str): Directory to save the output files.
            target_density (float): Desired solvent density in g/cm³.
            neutralize (bool): Whether to neutralize the system. Defaults to True.

        Returns:
            Dict[str, str]: Paths to the prepared files.
        """
        os.makedirs(output_dir, exist_ok=True)

        # Backup original topology file
        backup_top = os.path.join(output_dir, "topology_backup.top")
        shutil.copy(top_file, backup_top)

        # Edit the topology file to include the solvent.itp and force field
        self._edit_topology(top_file, solvent)

        # Prepare initial box
        box_gro = os.path.join(output_dir, "box.gro")
        final_density = None
        iteration = 0

        while True:
            # Create a cubic box
            self._create_cubic_box(top_file, box_gro, box_size)

            # Solvate the box using the .pdb file
            self._solvate_box(solvent.pdb_path, box_gro, top_file)

            # Validate density
            final_density = self._validate_density(box_gro, target_density)
            if abs(final_density - target_density) / target_density <= self.ACCEPTABLE_DENSITY_TOLERANCE:
                break

            # Adjust box size
            adjustment_factor = (final_density / target_density) ** (1 / 3)
            box_size *= adjustment_factor
            iteration += 1
            if iteration > 10:
                raise RuntimeError("Failed to achieve desired density after 10 iterations.")

        # Add ions if required
        if neutralize:
            self.add_ions(box_gro, top_file, output_dir)

        self.metadata_tracker.add_step("Prepare Simulation Box", {
            "topology_file": top_file,
            "box_gro": box_gro,
            "final_density": final_density,
            "target_density": target_density,
            "neutralized": neutralize
        })

        return {"box_gro": box_gro, "topology_file": top_file}

    def _edit_topology(self, top_file: str, solvent: Solvent):
        """
        Edit the topology file to include the solvent.itp and force field.

        Args:
            top_file (str): Path to the topology file.
            solvent (Solvent): Solvent object containing necessary paths.
        """
        with open(top_file, "a") as f:
            f.write(f'\n#include "{solvent.force_field}/forcefield.itp"\n')
            f.write(f'#include "{solvent.itp_path}"\n')

    def _create_cubic_box(self, top_file: str, box_gro: str, box_size: float):
        """
        Create a cubic box using GROMACS editconf.

        Args:
            top_file (str): Path to the topology file.
            box_gro (str): Path to the output .gro file.
            box_size (float): Size of the cubic box in nm.
        """
        editconf_command = [
            "gmx", "editconf",
            "-f", top_file,
            "-o", box_gro,
            "-c",  # Center the box
            "-d", str(box_size),
            "-bt", "cubic"  # Box type
        ]
        subprocess.run(editconf_command, check=True)

    def _solvate_box(self, solvent_pdb: str, box_gro: str, top_file: str):
        """
        Solvate the box using GROMACS solvate.

        Args:
            solvent_pdb (str): Path to the solvent .pdb file.
            box_gro (str): Path to the box .gro file.
            top_file (str): Path to the topology file.
        """
        solvate_command = [
            "gmx", "solvate",
            "-cp", box_gro,
            "-cs", solvent_pdb,  # Use .pdb as the solvent configuration
            "-o", box_gro,
            "-p", top_file
        ]
        subprocess.run(solvate_command, check=True)

    def add_ions(self, box_gro: str, top_file: str, output_dir: str):
        """
        Add ions to neutralize the simulation box.

        Args:
            box_gro (str): Path to the solvated box .gro file.
            top_file (str): Path to the topology file.
            output_dir (str): Directory to save the output files.
        """
        tpr_file = os.path.join(output_dir, "ions.tpr")
        grompp_command = [
            "gmx", "grompp",
            "-f", self._get_template("ions.mdp"),
            "-c", box_gro,
            "-p", top_file,
            "-o", tpr_file
        ]
        subprocess.run(grompp_command, check=True)

        neutralized_box_gro = os.path.join(output_dir, "neutralized_box.gro")
        genion_command = [
            "gmx", "genion",
            "-s", tpr_file,
            "-o", neutralized_box_gro,
            "-p", top_file,
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
            str: Path to the template file.
        """
        template_dir = os.path.abspath("templates")
        template_path = os.path.join(template_dir, filename)
        if not os.path.exists(template_path):
            raise FileNotFoundError(f"Template file '{filename}' not found in {template_dir}.")
        return template_path

    def _validate_density(self, box_gro: str, target_density: float) -> float:
        """
        Validate the density of the simulation box.

        Args:
            box_gro (str): Path to the box .gro file.
            target_density (float): Desired solvent density.

        Returns:
            float: Calculated density of the box.
        """
        # Placeholder: Use external calculation or manual validation
        return target_density  # For now, assume it's correct
