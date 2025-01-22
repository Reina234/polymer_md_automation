from pycgtool import CGSystem

cg = CGSystem("mapping_manual.txt")  # Load mapping file
cg.generate_cg_molecules(
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.tpr",
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.xtc",
)
