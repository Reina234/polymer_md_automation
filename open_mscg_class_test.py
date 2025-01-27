from open_mscg_topol_gen import OpenMSCGTopolGenerator
from open_mscg_pair import OpenMSCGForceMatcher

topol = OpenMSCGTopolGenerator("zzz/test_open_msg.yaml")
topol.create_topol("test.top")

fm = OpenMSCGForceMatcher(topol, "CG.lammpstrj")
fm.run("zzz/test")
