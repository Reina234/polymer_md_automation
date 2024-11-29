from preprocessing.template_utils import create_mdps, retrieve_mdps
from config.paths import TemplatedMdps

create_mdps(TemplatedMdps.NPT, 300)

print(retrieve_mdps(TemplatedMdps.NPT, 300))
