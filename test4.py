from preprocessing.file_converters.obabel_converter import OpenBabelConverter
from preprocessing.metadata_tracker import MetadataTracker
from preprocessing.parameterizers.acpype_parameterizer import ACPYPEParameterizer
from preprocessing.pdb_validators.base_pdb_validator import BasePDBValidator


metadata_tracker = MetadataTracker()
pdbtomol2converter = OpenBabelConverter(metadata_tracker)
print(pdbtomol2converter.input_file_type)

mol2_path = pdbtomol2converter.convert("styrene.pdb")

acpype = ACPYPEParameterizer(metadata_tracker)
acpype.parameterize(mol2_path, "test_new_struct")

validator = BasePDBValidator()
validator.skeletal_check("styrene.pdb")
