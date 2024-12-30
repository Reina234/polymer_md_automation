from A_modules.atomistic.gromacs.utils.utils import (
    get_gro_handler,
    calculate_minimum_box_size_from_df,
    export_gro_handler,
)

padding = 0.2
gro_file = "rdkit_test/POLY_GMX.gro"
gro_handler = get_gro_handler(gro_file)
box_size = calculate_minimum_box_size_from_df(gro_handler.content, padding=padding)
gro_handler.box_dimensions = box_size

export_gro_handler(gro_handler, output_path="trial.gro")
