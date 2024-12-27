from A_modules.atomistic.gromacs.utils.utils import (
    rename_residue_name_from_gro,
    prepare_solvent_topol,
    process_solvent_files,
)


# rename_residue_name_from_gro("ZZZZZ/test_output.gro", "TEST", "Z2")
# prepare_solvent_topol(
#    "TEST/POLY_GMX.top", "2", "TEST", "test/2/includes.itp", None, "Z2"
# )
process_solvent_files(
    "TEST/POLY_GMX.top", "TEST/POLY_GMX.itp", "TEST/POLY_GMX.gro", "TESTER", "Z2", "Z2"
)
