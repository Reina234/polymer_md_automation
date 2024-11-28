from gromacs.gromacs_utils import reformat_topol_file

topol_path = "output/test_new_struct/gromacs/topol.top"

test = reformat_topol_file(
    topol_path,
    "input/solvents/hexane.itp",
    "output/test_new_struct/acpype_output/POLY_GMX.itp",
    "input/forcefields/gromos54a7.ff",
    "test.top",
)
