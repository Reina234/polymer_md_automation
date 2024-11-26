from processing.pdb_validators.base_pdb_validator import BasePDBValidator

class PackmolPDBValidator(BasePDBValidator):
    """
    Validator for PDB files used with Packmol.
    """

    def validate(self) -> bool:
        """
        Validate the PDB file for Packmol compatibility.
        """
        # Perform the skeletal validation
        if not self._skeletal_check():
            return False

        # No additional Packmol-specific checks required
        return True

    @property
    def supports_fixing(self) -> bool:
        """
        Fixing is not required for Packmol.
        """
        return False
