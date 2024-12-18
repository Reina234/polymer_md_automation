from A_modules.shared.packmol.solvent_box import PackmolSolventBox
from data_models.solvent import Solvent
import os


def test_packmol_run():
    # Define test parameters
    solvent_file = "styrene.pdb"  # Path to a valid solvent PDB file
    output_dir = "test_output"
    output_file = "test_solvent_box.pdb"
    box_size_nm = [3, 3, 3]  # Example box size in nanometers
    solvent = Solvent("Hexane", 186.18, 660, solvent_file, "TMZK")

    # Create an instance of PackmolSolventBox
    packmol_operation = PackmolSolventBox()

    # Run the Packmol operation
    result = packmol_operation.run(
        solvent_file,
        output_file,
        output_dir="ZZZ",
        solvent=solvent,
        box_size_nm=box_size_nm,
    )

    # Validate the result
    assert os.path.exists(result), f"Output file not found: {result}"
    return result


result = test_packmol_run()
print(result)
from A_modules.shared.file_conversion.converters.editconf_pdb_to_gro import (
    EditconfPDBtoGROConverter,
)

converter = EditconfPDBtoGROConverter()
output_gro = converter.run(result, "ZZZ", "test_output.gro")
