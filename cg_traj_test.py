from mscg import Topology, Trajectory
import numpy as np

CG_traj = Trajectory("temp/cg_traj.lammpstrj", fmt="lammpstrj")
CG_traj.read_frame()
print(CG_traj.x)
print(CG_traj.x.shape)
print(CG_traj.box)

top = Topology.read_file("temp/cg_topol.top")
top.system_name = "cg_poly"

molecule_ids = []

for i in range(1, 1595):
    molecule_ids += [i, i, i, i, i, i]

top.save(
    "lammps",
    masses=[86.0, 105.0, 28.0, 104.0, 29.0],
    molecule=molecule_ids,
    box=np.vstack([np.zeros(3), CG_traj.box]).T,
    x=CG_traj.x,
)
