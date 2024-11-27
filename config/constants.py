from dataclasses import dataclass
from enum import Enum
AVOGADROS_NUMBER = 6.022e23
CM3_TO_NM3 = 1e21
NM3_TO_CM3 = 1e-21
ANGSTROM_TO_CM = 1e-8
DENSITY_TOLERANCE_PERCENTAGE = 5.0


class LengthUnits(Enum):
    ANGSTROM = "angstrom"
    NANOMETER = "nm"
    MICROMETER = "um"
    MILLIMETER = "mm"
    CENTIMETER = "cm"
    METER = "m"
    KILOMETER = "km"


CONVERSION_FACTORS_TO_M = {
    LengthUnits.ANGSTROM: 1e-10,
    LengthUnits.NANOMETER: 1e-9,
    LengthUnits.MICROMETER: 1e-6,
    LengthUnits.MILLIMETER: 1e-3,
    LengthUnits.CENTIMETER: 1e-2,
    LengthUnits.METER: 1,
    LengthUnits.KILOMETER: 1e3,
}

class MassUnits(Enum):
    AMU = "amu"
    GRAM = "g"
    KILOGRAM = "kg"

CONVERSION_FACTORS_TO_KG = {
    MassUnits.AMU: 1.66053906660e-27,
    MassUnits.GRAM: 1e-3,
    MassUnits.KILOGRAM: 1,
}