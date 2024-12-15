import numpy as np
from A_config.constants import (
    LengthUnits2,
    AVOGADROS_NUMBER,
    MassUnits2,
)
from typing import List


def calculate_minimum_box_size(
    atom_coordinates: List[List[float]],
    padding: float = 0.5 * 10**-9,
    units: LengthUnits2 = LengthUnits2.ANGSTROM,
    output_units: LengthUnits2 = LengthUnits2.METER,
) -> List[float]:

    conversion_factor = units.value
    coords = np.array(atom_coordinates)
    x_min, y_min, z_min = coords.min(axis=0) * conversion_factor
    x_max, y_max, z_max = coords.max(axis=0) * conversion_factor

    box_size = [
        (x_max - x_min) + 2 * padding,
        (y_max - y_min) + 2 * padding,
        (z_max - z_min) + 2 * padding,
    ]
    box_size = [dim / output_units.value for dim in box_size]
    return box_size


def calculate_density(
    molecular_weight: float,
    box_dimensions: List[float],
    length_units: LengthUnits2 = LengthUnits2.METER,
    mass_units: MassUnits2 = MassUnits2.GRAM,
    num_molecules: int = 1,
) -> float:
    """
    Calculate the density of the system based on box dimensions and molecular weight.

    Args:
        molecular_weight (float): Molecular weight of the molecule in g/mol.
        num_molecules (int): Number of molecules in the box.
        box_dimensions (List[float]): Box dimensions [x, y, z] in m.

    Returns:
        float: Density in kg/m³.
    """
    if not box_dimensions:
        raise ValueError("Box dimensions must be provided to calculate density.")

    volume_SI = np.prod(box_dimensions) * length_units.value**3
    molecular_weight_SI = molecular_weight * mass_units.value
    total_mass_SI = (num_molecules * molecular_weight_SI) / AVOGADROS_NUMBER

    density_SI = total_mass_SI / volume_SI
    return density_SI


def calculate_mass_per_unit_box(
    molecular_weight: float,
    mass_units: MassUnits2 = MassUnits2.GRAM,
    num_molecules: int = 1,
) -> float:
    molecular_weight_SI = molecular_weight * mass_units.value
    total_mass_SI = (num_molecules * molecular_weight_SI) / AVOGADROS_NUMBER
    return total_mass_SI


def calculate_volume_for_desired_density(
    molecular_weight: float,
    desired_density: float,
    mass_units: MassUnits2 = MassUnits2.GRAM,
    num_molecules: int = 1,
) -> float:
    """
    Calculate the required volume to achieve the desired density.

    Args:
        molecular_weight (float): Molecular weight of the molecule in g/mol.
        desired_density (float): Desired density in g/cm³.
        num_molecules (int): Number of molecules.

    Returns:
        float: Required volume in nm³.
    """

    total_mass_SI = calculate_mass_per_unit_box(
        molecular_weight, mass_units, num_molecules
    )

    volume_SI = total_mass_SI / desired_density
    return volume_SI


def scale_box_to_desired_volume(
    box_dimensions: List[float],
    desired_volume: float,
    output_box_units: LengthUnits2 = LengthUnits2.ANGSTROM,
    input_box_units: LengthUnits2 = LengthUnits2.METER,
    volume_units: LengthUnits2 = LengthUnits2.METER,
) -> List[float]:
    current_volume_SI = np.prod(box_dimensions) * input_box_units.value**3
    desired_volume_SI = desired_volume * volume_units.value**3
    scale_factor_SI = (desired_volume_SI / current_volume_SI) ** (1 / 3)
    if scale_factor_SI < 1:
        raise ValueError(
            f"Scale factor is less than 1 ({scale_factor_SI:.3f}). Desired volume "
            f"is smaller than the current box volume ({current_volume_SI:.3f} m³). "
            "This could lead to overlap or unrealistic density."
        )

    scaled_dimensions = [
        dim * scale_factor_SI / volume_units.value for dim in box_dimensions
    ]

    return scaled_dimensions


def calculate_num_particles(
    box_dimensions: List[float],
    molecular_weight: float,
    density_SI: float,
    box_units: LengthUnits2 = LengthUnits2.ANGSTROM,
    mass_units: MassUnits2 = MassUnits2.GRAM,
) -> int:

    volume_SI = np.prod(box_dimensions) * box_units.value**3
    molecular_weight_SI = molecular_weight * mass_units.value
    mass = volume_SI * density_SI
    num_particles = (mass / molecular_weight_SI) * AVOGADROS_NUMBER
    return round(int(num_particles))
