from collections import OrderedDict
from ZZZ_parser.section import Section
from ZZZ_parser.constants import GROMACS_CONSTRUCTS, GROMACSConstructs
from ZZZ_parser.handlers.base_handler import BaseHandler
from ZZZ_parser.handlers.default_handler import DefaultHandler
from typing import Dict, List, Optional, Tuple, OrderedDict


class FileSplitter:

    def __init__(
        self, full_constructs: Dict[str, GROMACSConstructs] = GROMACS_CONSTRUCTS
    ):
        self.full_constructs: Dict[str, GROMACSConstructs] = full_constructs
        self.suppressed_constructs: Optional[List] = None
        self._active_constructs: Dict[str, GROMACSConstructs] = full_constructs

    @property
    def active_constructs(self):
        # Dynamically filter suppressed patterns, this is to suppress checking for #includes and stuff within the if statements
        return {
            key: value
            for key, value in self.full_constructs.items()
            if not (self.suppressed_constructs and key in self.suppressed_constructs)
        }

    def split_file(self, filepath: str) -> OrderedDict[str, Section]:
        sections: OrderedDict[str, Section] = OrderedDict()
        current_section = Section(construct_type="default", handler_type=DefaultHandler)

        with open(filepath, "r") as file:
            for line in file:
                line = line.rstrip("\n")

                construct_type, name, handler = self._match_line(line)

                # Finalize the current section if a new construct type is matched
                if construct_type:
                    key = self._generate_key(
                        sections, current_section.construct_type, current_section.name
                    )
                    sections[key] = current_section

                    # Start a new section
                    current_section = Section(
                        construct_type=construct_type, handler_type=handler, name=name
                    )

                # Add line to the current section
                current_section.add_line(line)

        # Finalize the last section in the file
        if current_section.lines:
            key = self._generate_key(
                sections, current_section.construct_type, current_section.name
            )
            sections[key] = current_section

        return sections

    def _generate_key(
        self,
        sections: OrderedDict[str, Section],
        construct_type: str,
        name: Optional[str],
    ) -> str:
        """
        Generates a unique key for the section based on its type and name.
        """
        base_key = f"{construct_type}_{name or 'no_name'}"

        # Check for duplicates, but allow if it's the last active section
        if base_key in sections:
            last_key = next(reversed(sections), None)
            if last_key and last_key.startswith(base_key):
                return last_key  # Reuse the last key if it matches (active section)

            raise ValueError(f"Duplicate section found: {base_key}")

        return base_key

    def _match_line(self, line: str) -> Tuple[str, str, BaseHandler]:
        for construct_type, construct in self.active_constructs.items():
            match = construct.pattern.match(line)
            if match:
                name = match.group(1) if match.groups() else None
                self.suppressed_constructs = construct.suppress
                return construct_type, name, construct.handler

        return None, None, DefaultHandler
