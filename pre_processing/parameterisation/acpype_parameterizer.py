import os
import subprocess
from typing import Dict
from pre_processing.parameterisation.base_parameterizer import BaseParameterizer

class ACPYPEParameterizer(BaseParameterizer):
    """
    Handles parameterization using ACPYPE.
    """

    def __init__(self, metadata_tracker):
        super().__init__(metadata_tracker)

    def parameterize(self, input_file: str, output_dir: str) -> Dict[str, str]:
        """
        Parameterize the molecule using ACPYPE and update metadata.
        """
        # Ensure output directory exists
        parameterized_dir = os.path.join(output_dir, "parameterized")
        os.makedirs(parameterized_dir, exist_ok=True)

        # Construct ACPYPE command

    # Build the ACPYPE command
        command = [
            "acpype",
            "-i", input_file,       # Input file
            "-o", "gmx",            # Output format
            "-n", "0",              # Neutral charge
            "-a", "gaff2",           # Use GAFF2 force field
            "-b", "POLY"
        ]

        try:
            # Run the ACPYPE command
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"ACPYPE failed to parameterize {input_file}. Command: {' '.join(command)}. Error: {e}"
            )
        # Collect generated files
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        gro_file = os.path.join(parameterized_dir, f"{base_name}_GMX.gro")
        top_file = os.path.join(parameterized_dir, f"{base_name}_GMX.top")

        # Update metadata
        self.metadata_tracker.add_process_metadata(
            name="ACPYPE",
            version="2023.10",
            description="Generates GROMACS-compatible force field files.",
            options={
                "input_file": input_file,
                "output_dir": parameterized_dir,
                "gro_file": gro_file,
                "top_file": top_file,
                "charge": "neutral",
                "force_field": "GAFF"
            }
        )

        return {"gro_file": gro_file, "top_file": top_file}

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
                "force_field": "GAFF"
            }
        }
