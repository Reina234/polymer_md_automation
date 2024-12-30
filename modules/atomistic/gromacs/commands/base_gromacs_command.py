from abc import ABC, abstractmethod
from config.constants import LengthUnits
from modules.shared.command_line_operation import CommandLineOperation


class BaseGromacsCommand(CommandLineOperation, ABC):
    # output_name: str
    default_units: LengthUnits = LengthUnits.NANOMETER

    #    def __init_subclass__(cls):
    #        super().__init_subclass__()
    #        if not hasattr(cls, "output_name"):
    #            raise TypeError(f"Class {cls.__name__} must define 'output_name'.")

    @property
    def step_name(self) -> str:
        return "gromacs"

    @abstractmethod
    def run(self, *args, **kwargs):
        pass

    @abstractmethod
    def _create_command(self, *args, **kwargs):
        pass
