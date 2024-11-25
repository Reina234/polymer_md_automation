from pre_processing.converters.conversion_service import ConversionService
from pre_processing.converters.obabel_converter import OpenBabelConverter
from pre_processing.metadata.metadata_manager import MetadataManager
# Example: Preparing a PDB file from RDKit
from pdb_processing.providers.pdb_from_smiles_rdkit import PDBFromSmilesRDKit
from pdb_processing.managers.pdb_manager import PDBManager
from solvent.solvent import Solvent
from solvent.packmol_solvation import PackmolBoxGenerator
from force_field_handler import   ForceFieldHandler
from pre_processing.parameterisation.acpype_parameterizer import ACPYPEParameterizer
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

metadata_manager = MetadataManager()

conversion_service = ConversionService()
conversion_service.register_converter(OpenBabelConverter("pdb", "mol2", metadata_manager))
parameterizer = ACPYPEParameterizer(conversion_service, metadata_manager)
parameterization_dir = parameterizer.parameterize(pdb_path, output_dir)

    # Add user notes
metadata_manager.add_notes("This workflow was executed for testing purposes.")

    # Save metadata
metadata_manager.save_metadata(output_dir)
