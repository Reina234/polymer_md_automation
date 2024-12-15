from abc import ABC, abstractmethod
from A_modules.shared.command_line_operation import CommandLineOperation
from A_modules.shared.metadata_tracker import MetadataTracker
from typing import Optional


class BaseConverter(CommandLineOperation, ABC):
    input_file_type: str
    output_file_type: str
    program: str

    @property
    def step_name(self) -> str:
        return "conversion"

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if not hasattr(cls, "input_file_type") or not hasattr(cls, "output_file_type"):
            raise TypeError(
                f"Class {cls.__name__} must define 'input_file_type' and 'output_file_type'."
            )

    def __init__(self, metadata_tracker: Optional[MetadataTracker] = None):
        super().__init__(metadata_tracker)

    @abstractmethod
    def run(
        self,
        input_file_path: str,
        output_dir: Optional[str] = None,
        additional_notes: Optional[str] = None,
        verbose: bool = False,
    ) -> str:
        pass

    def metadata(
        self,
        input_file_path: str,
        output_file_path: str,
        additional_notes: Optional[str] = None,
    ) -> dict:
        """
        Generate metadata for the conversion step.

        :param input_file_path: input file path
        :type input_file_path: str
        :param output_file_path: output file path
        :type output_file_path: str
        :param additional_notes: any additional notes defaults to None
        :type additional_notes: Optional[str], optional
        :return: metadata
        :rtype: Dict
        """
        return {
            "program(s) used": "program",
            "details": f"converts {input_file_path} to {output_file_path}",
            "action(s)": f"converted file type {self.input_file_type} to {self.output_file_type}",
            "additional notes": additional_notes,
        }
