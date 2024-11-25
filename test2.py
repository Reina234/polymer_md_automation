from processing.atomistic.parameterizers.acpype_parameterizer import ACPYPEParameterizer
from processing.common.metadata_tracker import MetadataTracker
tracker = MetadataTracker()
acpype = ACPYPEParameterizer(tracker)
pdb_path = "hexane.pdb"

acpype.parameterize(pdb_path, "test_output")