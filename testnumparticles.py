from A_modules.shared.utils.calculation_utils import calculate_num_particles
from A_config.constants import LengthUnits2

print(calculate_num_particles([10, 10, 10], 102, 990, LengthUnits2.NANOMETER))
