import os
import subprocess
import shutil
import logging
from typing import Dict
from processing.atomistic.parameterizers.base_parameterizer import BaseParameterizer
from processing.metadata_tracker import MetadataTracker
from processing.file_converters.obabel_converter import OpenBabelConverter

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


class ACPYPEParameterizer(BaseParameterizer):
    """
    Parameterizes a polymer using ACPYPE for atomistic simulations.
    """

    BASE_TEMPORARY_DIRECTORY = "temp/atomistic"
    BASE_SUBDIR = "acpype_gmx_files"

    def __init__(self, metadata_tracker: MetadataTracker):
        super().__init__(metadata_tracker)
        self.converter = OpenBabelConverter("pdb", "mol2", metadata_tracker)

    def parameterize(self, input_file: str, output_dir: str) -> Dict[str, str]:
        """
        Parameterize the polymer using ACPYPE.

        Args:
            input_file (str): Full path to the input .pdb file.
            output_dir (str): Directory to save parameterized output files.

        Returns:
            Dict[str, str]: Paths to the generated .gro, .top, and .itp files.
        """
        logger.info("Starting ACPYPE parameterization.")
        logger.debug(f"Input file: {input_file}")
        logger.debug(f"Output directory: {output_dir}")

        if not input_file.endswith(".pdb"):
            raise ValueError(f"Input file must be a .pdb file: {input_file}")
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")

        # Ensure directories are absolute
        input_file = os.path.abspath(input_file)
        output_dir = os.path.abspath(output_dir)
        os.makedirs(output_dir, exist_ok=True)

        polymer_name = os.path.splitext(os.path.basename(input_file))[0]
        temp_dir = os.path.abspath(self.BASE_TEMPORARY_DIRECTORY)
        #temp_dir = os.path.abspath(os.path.join(self.BASE_TEMPORARY_DIRECTORY, polymer_name))
        os.makedirs(temp_dir, exist_ok=True)

        # Step 1: Convert PDB to MOL2 using OpenBabelConverter
        mol2_file = os.path.join(output_dir, f"{polymer_name}.mol2")
        if not os.path.exists(mol2_file):
            mol2_file = self.converter.convert(input_file, output_dir)

        logger.info(f"Converted PDB to MOL2: {mol2_file}")
        if not os.path.exists(mol2_file):
            raise FileNotFoundError(f"Conversion to MOL2 failed: {mol2_file}")

        # Step 2: Run ACPYPE
        base_name = "POLY"
        acpype_command = [
            "acpype",
            "-i", mol2_file,
            "-o", "gmx",
            "-n", "0",
            "-a", "gaff2",
            "-b", base_name
        ]
        acpype_temp_dir = os.path.join(temp_dir, f"{base_name}.acpype")
        os.makedirs(acpype_temp_dir, exist_ok=True)
        logger.info("Running ACPYPE...")
        logger.debug(f"ACPYPE command: {' '.join(acpype_command)}")

        try:
            result = subprocess.run(
                acpype_command,
                cwd=temp_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logger.debug(f"ACPYPE stdout: {result.stdout}")
            logger.debug(f"ACPYPE stderr: {result.stderr}")

            if result.returncode != 0:
                raise RuntimeError(f"ACPYPE failed:\nSTDOUT: {result.stdout}\nSTDERR: {result.stderr}")

            logger.info("ACPYPE completed successfully.")
            logger.debug(f"ACPYPE Temporary Directory Contents: {os.listdir(acpype_temp_dir)}")

            # Step 3: Locate output files
            gro_file = os.path.join(acpype_temp_dir, f"{base_name}_GMX.gro")
            top_file = os.path.join(acpype_temp_dir, f"{base_name}_GMX.top")
            itp_file = os.path.join(acpype_temp_dir, f"{base_name}_GMX.itp")

            if not all(os.path.exists(file) for file in [gro_file, top_file, itp_file]):
                raise FileNotFoundError(f"One or more ACPYPE output files are missing in {acpype_temp_dir}")

            output_gro = os.path.join(output_dir, self.BASE_SUBDIR, f"{polymer_name}.gro")
            output_top = os.path.join(output_dir, self.BASE_SUBDIR, f"{polymer_name}.top")
            output_itp = os.path.join(output_dir, self.BASE_SUBDIR, f"{polymer_name}.itp")

            os.makedirs(os.path.dirname(output_gro), exist_ok=True)

            for src, dest in zip([gro_file, top_file, itp_file], [output_gro, output_top, output_itp]):
                logger.info(f"Moving {src} to {dest}.")
                shutil.move(src, dest)

            # Clean up temporary files
            logger.info("Cleaning up temporary files.")
            shutil.rmtree(acpype_temp_dir)

            # Update metadata
            logger.info("Updating metadata for ACPYPE parameterization.")
            self.metadata_tracker.add_step("ACPYPE Parameterization", {
                "input_file": input_file,
                "output_files": [output_gro, output_top, output_itp],
                "mol2_file": mol2_file,
                "command": ' '.join(acpype_command),
                "force_field": "GAFF2"
            })

            logger.info("ACPYPE parameterization completed successfully.")
            return {
                "gro_file": output_gro,
                "top_file": output_top,
                "itp_file": output_itp,
                "mol2_file": mol2_file
            }

        except Exception as e:
            logger.exception("Error occurred during ACPYPE parameterization.")
            raise RuntimeError(f"ACPYPE failed to parameterize polymer. Error: {str(e)}")

    def metadata(self) -> Dict:
        """
        Return metadata about ACPYPE.
        """
        return {
            "name": "ACPYPE Parameterizer",
            "version": "2023.10",
            "description": "Generates GROMACS-compatible force field files.",
            "options": {
                "charge": "neutral",
                "force_field": "GAFF2"
            }
        }
