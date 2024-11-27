from processing.atomistic.parameterizers.acpype_parameterizer import ACPYPEParameterizer
from processing.metadata_tracker import MetadataTracker
tracker = MetadataTracker()
acpype = ACPYPEParameterizer(tracker)
pdb_path = "styrene.pdb"

acpype.parameterize(pdb_path, "test_output_styrene")