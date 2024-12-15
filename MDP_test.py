from A_modules.atomistic.gromacs.utils.mdp_utils import create_mdps
from A_modules.atomistic.gromacs.gromacs_config import TemplatedMdps

create_mdps(TemplatedMdps.NVT, "TEST", 300)
