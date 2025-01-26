from mscg import *
from mscg.cli import cgmap

cgmap.main(map="zzz/test_open_msg.yaml", traj="temp/production.trr", out="CG.lammpstrj")

CG_traj = Trajectory("CG.lammpstrj", fmt="lammpstrj")
CG_traj.read_frame()
CG_traj.x.shape
