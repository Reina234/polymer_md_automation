from data_models.solvent import Solvent
from processing.metadata_tracker import MetadataTracker
from processing.pdb_validators.gromacs_pdb_validator import GROMACSPDBValidator
from processing.parsers.EDIT_pdb_parser import PDBParser

solvent = Solvent(
    "Hexane",
    86.18,
    0.660,
    "input/solvents/itp/hexane.itp",
    "input/solvents/pdb/hexane.pdb",
    "gromos54a7.ff",
)
metadata_tracker = MetadataTracker()
pdb_validator = GROMACSPDBValidator()


# path = pdb_validator.validate("test_pdb.pdb", solvent, replace_box = True)
#
# print(path)

pdb_parser = PDBParser("test_pdb.pdb")


solvent.pdb_path = "test_pdb.pdb"
solvent_gro_path = "test_output/acpype_gmx_files/hexane.gro"

polymer_gro_file = "test_output_styrene/acpype_gmx_files/styrene.gro"
polymer_top_file = "test_output_styrene/acpype_gmx_files/styrene.top"
output_dir = "test_output_styrene"
from gromacs.solvation.box_creator import BoxCreator
from gromacs.solvation.NOTE_solvate_box import SolvateBox
from gromacs.solvation.ion_adder import IonAdder

# Initialize the box preparer
box_preparer = BoxCreator(metadata_tracker)

# NOTE: need to add metadata tracker?
box_size = 3

print("[INFO] Creating a cubic simulation box...")
box_file = box_preparer.create_box(
    input_gro_path=polymer_gro_file, output_dir=output_dir, box_size=box_size
)
print(f"[INFO] Cubic box created: {box_file}")


print("[INFO] Solvating the cubic box...")
solvate_box = SolvateBox(metadata_tracker)

solvated_box_results = solvate_box.solvate(
    box_file=box_file,
    solvent_file=solvent_gro_path,
    polymer_top_file=polymer_top_file,
    output_dir=output_dir,
)

print(f"[INFO] Solvated box created: {solvated_box_results['solvated_gro_file']}")

print("[INFO] Adding ions to the solvated box...")
ion_adder = IonAdder(metadata_tracker)
neutralized_results = ion_adder.add_ions(
    input_gro=solvated_box_results["solvated_gro_file"],
    topology_file=solvated_box_results["topology_file"],
    output_dir=output_dir,
)

print(f"[INFO] Ion addition completed. Results: {neutralized_results}")
