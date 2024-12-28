from A_modules.atomistic.gromacs.utils.utils import (
    rename_residue_name_from_gro,
    prepare_solvent_topol,
    process_solvent_files,
    process_solvent_itp,
)


# rename_residue_name_from_gro("ZZZZZ/test_output.gro", "TEST", "Z2")
# prepare_solvent_topol(
#    "TEST/POLY_GMX.top", "2", "TEST", "test/2/includes.itp", None, "Z2"
# )

itp_file = process_solvent_itp("TEST_RUN_27/POLY_GMX.itp", "TMZK", "Z2")
process_solvent_files(
    "TEST_RUN_27/POLY_GMX.top",
    itp_file,
    "TEST_RUN_27/POLY_GMX.gro",
    new_residue_name="TEST",
    output_gro_dir="Z2",
    output_topol_dir="Z2",
)
