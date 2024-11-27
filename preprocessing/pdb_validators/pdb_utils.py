import numpy as np
from config.constants import (
    LengthUnits,
    AVOGADROS_NUMBER,
    CONVERSION_FACTORS_TO_M,
    MassUnits,
    CONVERSION_FACTORS_TO_KG,
)
from typing import List


def calculate_minimum_box_size(
    atom_coordinates: List[List[float]],
    padding: float = 10**-9,
    units: LengthUnits = LengthUnits.ANGSTROM,
) -> List[float]:
    """
    Calculate the minimum bounding box size based on atom positions.

    Args:
        atom_coordinates (List[List[float]]): List of atom coordinates [x, y, z].
        padding (float): Padding to add to each dimension - padding is in m
        units = units of file !!! need to edit this

    Returns:
        List[float]: Minimum box dimensions [x, y, z] in nm.
    """
    conversion_factor = CONVERSION_FACTORS_TO_M[units]
    coords = np.array(atom_coordinates)
    x_min, y_min, z_min = coords.min(axis=0) * conversion_factor
    x_max, y_max, z_max = coords.max(axis=0) * conversion_factor

    box_size = [
        (x_max - x_min) + 2 * padding,
        (y_max - y_min) + 2 * padding,
        (z_max - z_min) + 2 * padding,
    ]
    return box_size


def calculate_density(
    molecular_weight: float,
    box_dimensions: List[float],
    length_units: LengthUnits = LengthUnits.METER,
    mass_units: MassUnits = MassUnits.GRAM,
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

    volume_SI = np.prod(box_dimensions) * CONVERSION_FACTORS_TO_M[length_units] ** 3
    molecular_weight_SI = molecular_weight * CONVERSION_FACTORS_TO_KG[mass_units]
    total_mass_SI = (num_molecules * molecular_weight_SI) / AVOGADROS_NUMBER

    density_SI = total_mass_SI / volume_SI
    return density_SI


def calculate_mass_per_unit_box(
    molecular_weight: float,
    mass_units: MassUnits = MassUnits.GRAM,
    num_molecules: int = 1,
) -> float:
    molecular_weight_SI = molecular_weight * CONVERSION_FACTORS_TO_KG[mass_units]
    total_mass_SI = (num_molecules * molecular_weight_SI) / AVOGADROS_NUMBER
    return total_mass_SI


def calculate_volume_for_desired_density(
    molecular_weight: float,
    desired_density: float,
    mass_units: MassUnits = MassUnits.GRAM,
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
    output_box_units: LengthUnits = LengthUnits.ANGSTROM,
    input_box_units: LengthUnits = LengthUnits.METER,
    volume_units: LengthUnits = LengthUnits.METER,
) -> List[float]:
    current_volume_SI = (
        np.prod(box_dimensions) * CONVERSION_FACTORS_TO_M[input_box_units] ** 3
    )
    desired_volume_SI = desired_volume * CONVERSION_FACTORS_TO_M[volume_units] ** 3
    scale_factor_SI = (desired_volume_SI / current_volume_SI) ** (1 / 3)
    if scale_factor_SI < 1:
        raise ValueError(
            f"Scale factor is less than 1 ({scale_factor_SI:.3f}). Desired volume "
            f"is smaller than the current box volume ({current_volume_SI:.3f} m³). "
            "This could lead to overlap or unrealistic density."
        )

    scaled_dimensions = [
        dim * scale_factor_SI / CONVERSION_FACTORS_TO_M[output_box_units]
        for dim in box_dimensions
    ]

    return scaled_dimensions
