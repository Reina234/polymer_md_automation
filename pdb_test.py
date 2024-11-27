# Example: Preparing a PDB file from RDKit
from pdb_processing.providers.pdb_from_smiles_rdkit import PDBFromSmilesRDKit
from pdb_processing.managers.pdb_manager import PDBManager
from data_models.solvent import Solvent
from solvent.packmol_solvation import PackmolBoxGenerator
from force_field_handler import   ForceFieldHandler
# Define solvent
ethanol_solvent = Solvent(name="Ethanol", density=0.789, molecular_weight=46.07)

# Create provider
rdkit_provider = PDBFromSmilesRDKit(smiles="CCO", solvent = ethanol_solvent, additional_notes="Generated for ethanol test suite.")

# Initialize manager
pdb_manager = PDBManager(provider=rdkit_provider)

# Prepare PDB
output_dir = "output/pdb_files"
pdb_path = pdb_manager.prepare_pdb(output_dir,  identifier="Test")

# Retrieve metadata
metadata = pdb_manager.metadata
print(f"Metadata:\n{metadata}")

# Retrieve PDB path
print(f"PDB Path: {pdb_manager.pdb_path}")


#packmol_gen = PackmolBoxGenerator(pdb_manager=pdb_manager, box_size=30, output_dir="output/packmol_files")
#packmol_gen.generate_box()

forcefieldhandler = ForceFieldHandler()
forcefieldhandler.parameterize_with_acpype(pdb_manager.pdb_path, "output")
