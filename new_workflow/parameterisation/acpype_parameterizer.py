class ACPYPEParameterizer(BaseParameterizer):
    """
    Parameterizes molecules using ACPYPE with the GAFF force field.
    """

    def __init__(self, file_converter: FileConverter):
        super().__init__("ACPYPE", file_converter)

    def parameterize(self, input_file: str, output_dir: str) -> str:
        """
        Parameterize the molecule using ACPYPE.
        Converts to MOL2 if needed and generates GROMACS-compatible files.
        """
        # Convert to MOL2 if input is not already MOL2
        if not input_file.endswith(".mol2"):
            input_file = self.file_converter.convert(input_file, "mol2", output_dir)

        acpype_dir = os.path.join(output_dir, os.path.splitext(os.path.basename(input_file))[0] + "_acpype")
        subprocess.run(
            [
                "acpype",
                "-i", input_file,
                "-o", "gmx",  # GROMACS-compatible output
                "-n", "0"     # Neutral charge
            ],
            check=True
        )
        return acpype_dir
