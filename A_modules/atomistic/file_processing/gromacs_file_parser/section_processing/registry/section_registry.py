from typing import List, Dict
from A_modules.atomistic.file_processing.gromacs_file_parser.section_processing.registry.prebuilt_registries import (
    acpype_generated_sections,
)


class SectionRegistry:
    def __init__(self, prebuilt_registry: Dict[str, List] = acpype_generated_sections):
        self.registry = prebuilt_registry

    def register_section(self, section_name: str, headers: List[str]):
        """
        Registers a section and its expected headers.
        """
        if section_name in self.registry:
            raise ValueError(f"Section '{section_name}' is already registered.")
        self.registry[section_name] = headers
        # is this needed??? NOTE: probably not?

    def get_headers(self, section_name: str) -> List[str]:
        """
        Retrieves the expected headers for a given section.
        """
        if section_name not in self.registry:
            raise ValueError(f"Section '{section_name}' is not registered.")
        return self.registry[section_name]

    def validate_headers(self, section_name: str, provided_headers: List[str]):
        """
        Validates that the provided headers match the registered headers.
        """
        expected_headers = self.get_headers(section_name)
        if expected_headers != provided_headers:
            raise ValueError(
                f"Headers for section '{section_name}' do not match. "
                f"Expected: {expected_headers}, Provided: {provided_headers}"
            )
