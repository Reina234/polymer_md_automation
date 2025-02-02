from abc import ABC, abstractmethod
from typing import Optional, Dict
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BaseMoltemplateMolecule(ABC):
    def __init__(self, open_mscg_mass_map: Optional[Dict[str, float]] = None):
        self.premade_mass_map = open_mscg_mass_map
        self.mass_mapping = {}

    def _compare_mass_maps(self) -> None:
        mass_map_1 = self.mass_mapping
        mass_map_2 = self.premade_mass_map
        if not mass_map_2 or not mass_map_1:
            logger.warning(
                "No mass map provided for comparison. Skipping mass map comparison."
            )
            return
        common_beads = set(mass_map_1.keys()) & set(
            mass_map_2.keys()
        )  # Find shared keys

        for bead_type in common_beads:
            mass1, mass2 = mass_map_1[bead_type], mass_map_2[bead_type]

            # Check if values match when rounded to the nearest integer
            if round(mass1) != round(mass2):
                logger.warning(
                    f"Masses for bead type {bead_type} do not match. For mass_map_1: {mass1}, for mass_map_2: {mass2}"
                )

    @abstractmethod
    def generate_lt(self) -> str:
        pass
