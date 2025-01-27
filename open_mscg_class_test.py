from modules.open_mscg.topol_generator import OpenMSCGTopolGenerator
from modules.open_mscg.force_matcher import OpenMSCGForceMatcher

topol = OpenMSCGTopolGenerator("zzz/test_open_msg.yaml")
topol.create_topol("test.top")

fm = OpenMSCGForceMatcher(topol, "CG.lammpstrj")
fm.run_cgfm("zzz/test")
fm.run_cgdump(output_dir="zzz")
fm.visualise_tables("zzz")
