import MDAnalysis as mda
from hoomd import data
import gsd.hoomd

# Load the GROMACS atomistic trajectory
u = mda.Universe(
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.gro",
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.xtc",
)

# Create GSD file
with gsd.hoomd.open("cg_trajectory.gsd", "w") as gsd_file:
    for ts in u.trajectory:
        frame = data.frame_from_universe(u)
        gsd_file.append(frame)
