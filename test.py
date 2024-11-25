from processing.common.file_converters.REDUNDANT_conversion_service import ConversionService
from processing.common.file_converters.obabel_converter import OpenBabelConverter
from processing.common.metadata_tracker import MetadataTracker
# Example: Preparing a PDB file from RDKit
from pdb_processing.providers.pdb_from_smiles_rdkit import PDBFromSmilesRDKit
from pdb_processing.managers.pdb_manager import PDBManager
from data_models.solvent import Solvent
from solvent.packmol_solvation import PackmolBoxGenerator

from processing.atomistic.parameterizers.acpype_parameterizer import ACPYPEParameterizer
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


# Initialize MetadataTracker
metadata_tracker = MetadataTracker()

# Setup ConversionService and register converters
conversion_service = ConversionService()
conversion_service.register_converter(OpenBabelConverter("pdb", "mol2", metadata_tracker))
conversion_service.register_converter(OpenBabelConverter("mol2", "pdb", metadata_tracker))

converter = OpenBabelConverter("pdb", "mol2", metadata_tracker)

mol2_file = converter.convert(pdb_path, output_dir)

    # Step 2: Parameterize using ACPYPE
parameterizer = ACPYPEParameterizer(metadata_tracker)
parameterized_files = parameterizer.parameterize(mol2_file, output_dir, "ethanol")

    # Save final metadata
output_dir = "output"
metadata_tracker.save_metadata(output_dir)

print(f"Parameterization complete. Files saved: {parameterized_files}")
print(f"Metadata saved to: {output_dir}/metadata.json")
