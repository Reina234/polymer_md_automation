from A_modules.shared.utils.calculation_utils import calculate_num_particles

num_molecules = calculate_num_particles(
    box_dimensions=[30, 30, 30],
    molecular_weight=100,
    density_SI=660,
)
print(num_molecules)
