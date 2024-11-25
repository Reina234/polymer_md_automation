from pre_processing.converters.conversion_service import ConversionService
from pre_processing.metadata.metadata_manager import MetadataManager
from pre_processing.parameterisation.base_parameterizer import BaseParameterizer
import os
import subprocess

class ACPYPEParameterizer(BaseParameterizer):
    """
    Parameterizes molecules using ACPYPE with the GAFF force field.
    """

    def __init__(self, conversion_service: ConversionService, metadata_manager: MetadataManager):
        super().__init__("ACPYPE", conversion_service, metadata_manager)

    def parameterize(self, input_file: str, output_dir: str) -> str:
        input_file = self.conversion_service.convert(input_file, "mol2", output_dir)
        acpype_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(input_file))[0] + "_acpype")
        subprocess.run(
            [
                "acpype",
                "-i", input_file,
                "-b", "POLY",
                "-o", "gmx",
                "-n", "0"
            ],
            check=True
        )
        return acpype_dir

    def metadata(self) -> dict:
        return {
            "name": "ACPYPE Parameterizer",
            "version": "2023.1",
            "description": "Generates GROMACS-compatible force field files using the GAFF force field.",
            "options": {
                "charge": "Neutral (0)",
                "force_field": "GAFF"
            }
        }
