from data_models.solvent import Solvent
from processing.metadata_tracker import MetadataTracker
from processing.pdb_validators.gromacs_pdb_validator import GROMACSPDBValidator 

solvent = Solvent("Hexane", 86.18, 0.660, "input/solvents/itp/hexane.itp", "input/solvents/pdb/hexane.pdb", "gromos54a7.ff")
metadata_tracker = MetadataTracker()
pdb_validator = GROMACSPDBValidator(solvent = solvent)

path = pdb_validator.validate("test_pdb.pdb")
print(path)




solvent.pdb_path = "test_pdb.pdb"


polymer_gro_file = "test_output_styrene/acpype_gmx_files/styrene.gro"
polymer_top_file = "test_output_styrene/acpype_gmx_files/styrene.top"
output_dir = "test_output_styrene"
from gromacs.solvation.box_creator import BoxCreator
from gromacs.solvation.NOTE_solvate_box import SolvateBox
from gromacs.solvation.ion_adder import IonAdder

# Initialize the box preparer
box_preparer = BoxCreator(metadata_tracker)

#NOTE: need to add metadata tracker?
box_size = 3

print("[INFO] Creating a cubic simulation box...")
box_file = box_preparer.create_box(input_gro = polymer_gro_file, output_dir = output_dir, box_size = box_size)
print(f"[INFO] Cubic box created: {box_file}")


print("[INFO] Solvating the cubic box...")
solvate_box = SolvateBox(metadata_tracker)

solvated_box_results = solvate_box.solvate(box_file = box_file,
                                           solvent = solvent, 
                                           polymer_top_file = polymer_top_file,
                                           output_dir = output_dir)

print(f"[INFO] Solvated box created: {solvated_box_results['solvated_gro']}")