from preprocessing.file_converters.obabel_converter import OpenBabelConverter
from preprocessing.metadata_tracker import MetadataTracker
from preprocessing.parameterizers.acpype_parameterizer import ACPYPEParameterizer
from preprocessing.validators.pdb_validators.base_pdb_validator import BasePDBValidator
from preprocessing.validators.pdb_validators.density_gromacs_pdb_validator import (
    GROMACSPDBValidator,
)

metadata_tracker = MetadataTracker()
pdbtomol2converter = OpenBabelConverter(metadata_tracker)
print(pdbtomol2converter.input_file_type)

mol2_path = pdbtomol2converter.convert("styrene.pdb")

acpype = ACPYPEParameterizer(metadata_tracker)
acpype.parameterize(mol2_path, "test_new_struct")

validator = GROMACSPDBValidator(MetadataTracker())
validator.validate("styrene.pdb", solvent, output_file_path="styrene.pdb")
