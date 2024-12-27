from A_modules.atomistic.gromacs.utils.utils import (
    rename_residue_name_from_gro,
    modify_topol_file,
)


# rename_residue_name_from_gro("ZZZZZ/test_output.gro", "TEST", "Z2")
modify_topol_file("TEST/POLY_GMX.top", "TEST", "2", "test/2/includes.itp", None, "Z2")
