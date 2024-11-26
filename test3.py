from data_models.solvent import Solvent
from processing.metadata_tracker import MetadataTracker
from simulation_scripts.atomistic.solvation_box_preparer import SolvationBoxPreparer
from processing.pdb_validators.gromacs_pdb_validator import GROMACSPDBValidator 

solvent = Solvent("Hexane", 86.18, 0.660, "input/solvents/itp/hexane.itp", "input/solvents/pdb/hexane.pdb", "gromos54a7.ff")
metadata_tracker = MetadataTracker()
pdb_validator = GROMACSPDBValidator(solvent, metadata_tracker)
pdb_validator.validate("test_pdb.pdb")




solvent.pdb_path = "test_pdb.pdb"

# Initialize the box preparer
box_preparer = SolvationBoxPreparer(metadata_tracker)

# Prepare the simulation box
results = box_preparer.prepare_box(
    solvent=solvent,
    polymer_top_file="test_output/acpype_gmx_files/hexane.top",
    polymer_gro_file="test_output/acpype_gmx_files/hexane.gro",
    box_size=5.0,  # nm
    output_dir="output/simulation_box",
    neutralize=True
)
