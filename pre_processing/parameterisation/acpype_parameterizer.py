import os
import shutil
import subprocess
from typing import Dict
from pre_processing.parameterisation.base_parameterizer import BaseParameterizer


class ACPYPEParameterizer(BaseParameterizer):
    """
    Handles parameterization using ACPYPE.
    """

    def __init__(self, metadata_tracker):
        super().__init__(metadata_tracker)

    def parameterize(self, input_file: str, output_dir: str, polymername: str) -> Dict[str, str]:
        """
        Parameterize the molecule using ACPYPE, save key files, and clean up.
        Args:
            input_file (str): Path to the input file (e.g., MOL2).
            output_dir (str): Directory to store parameterized output files.
            polymername (str): Name of the polymer for file naming.
        Returns:
            Dict[str, str]: Paths to the generated .itp, .top, and .gro files.
        """
        # Ensure absolute paths
        input_file = os.path.abspath(input_file)
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        # Define ACPYPE base name
        base_name = "POLY"

        try:
            # Build the ACPYPE command
            command = [
                "acpype",
                "-i", input_file,  # Input file
                "-o", "gmx",       # Output format
                "-n", "0",         # Neutral charge
                "-a", "gaff2",     # Use GAFF2 force field
                "-b", base_name    # Base name for output files
            ]

            # Run the ACPYPE command directly in the current working directory
            result = subprocess.run(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            if result.returncode != 0:
                # Print the stderr for debugging
                print("ACPYPE stderr:", result.stderr)
                print("ACPYPE stdout:", result.stdout)
                raise RuntimeError(f"ACPYPE failed with error: {result.stderr}")

            # Locate the ACPYPE output folder
            acpype_folder = f"{base_name}.acpype"

            if not os.path.exists(acpype_folder):
                raise FileNotFoundError(f"ACPYPE output folder '{acpype_folder}' not found.")

            # Paths to the expected files in the ACPYPE folder
            gro_file = os.path.join(acpype_folder, f"{base_name}_GMX.gro")
            top_file = os.path.join(acpype_folder, f"{base_name}_GMX.top")
            itp_file = os.path.join(acpype_folder, f"{base_name}_GMX.itp")

            # Validate file existence
            if not (os.path.exists(gro_file) and os.path.exists(top_file) and os.path.exists(itp_file)):
                raise FileNotFoundError("Expected ACPYPE output files were not generated in the ACPYPE folder.")

            # Rename and move extracted files to the output directory
            renamed_gro = os.path.join(output_dir, f"{polymername}_truncated.gro")
            renamed_top = os.path.join(output_dir, f"{polymername}_truncated.top")
            renamed_itp = os.path.join(output_dir, f"{polymername}_truncated.itp")

            shutil.move(gro_file, renamed_gro)
            shutil.move(top_file, renamed_top)
            shutil.move(itp_file, renamed_itp)

            # Update metadata
            self.metadata_tracker.add_process_metadata(
                name="ACPYPE",
                version="2023.10",
                description="Generates GROMACS-compatible force field files.",
                options={
                    "input_file": input_file,
                    "output_dir": output_dir,
                    "gro_file": renamed_gro,
                    "top_file": renamed_top,
                    "itp_file": renamed_itp,
                    "charge": "neutral",
                    "force_field": "GAFF2"
                }
            )

            return {
                "gro_file": renamed_gro,
                "top_file": renamed_top,
                "itp_file": renamed_itp
            }

        finally:
            # Clean up the ACPYPE folder and any remaining temporary files
            acpype_folder = f"{base_name}.acpype"
            if os.path.exists(acpype_folder):
                shutil.rmtree(acpype_folder)

            for file in os.listdir():
                if file.startswith(base_name):
                    file_path = os.path.join(os.getcwd(), file)
                    if os.path.isfile(file_path):
                        os.remove(file_path)

    def metadata(self) -> Dict:
        """
        Return metadata for the ACPYPE parameterizer.
        """
        return {
            "name": "ACPYPE Parameterizer",
            "version": "2023.10",
            "description": "Generates GROMACS-compatible force field files using the GAFF force field.",
            "options": {
                "charge": "neutral",
                "force_field": "GAFF2"
            }
        }
