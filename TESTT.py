from A_modules.atomistic.gromacs.workflows.AA_packmol import (
    create_density_accurate_solvent_box,
)

solvent_file = "styrene.pdb"
output_file = "solvent_box.gro"
box_size = [50.0, 50.0, 50.0]  # Box dimensions in Å
target_density = 0.3  # g/cm³ for water
molecular_weight = 80.015  # g/mol for water

try:
    solvent_box = create_density_accurate_solvent_box(
        solvent_file=solvent_file,
        output_file=output_file,
        box_size=box_size,
        target_density=target_density,
        molecular_weight=molecular_weight,
    )
    print(f"Solvent box created: {solvent_box}")
except RuntimeError as e:
    print(f"Error: {e}")
