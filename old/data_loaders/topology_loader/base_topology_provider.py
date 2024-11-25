# core/topology_provider.py
from abc import ABC, abstractmethod
from typing import Dict


class TopologyProvider(ABC):
    """
    Abstract base class for providing topology and force field parameters.
    """

    @abstractmethod
    def get_topology(self) -> Dict:
        """
        Generate or load topology/parameters and return as a dictionary.
        :return: A dictionary containing topology and parameters.
        """
        pass
