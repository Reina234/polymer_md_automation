from collections import OrderedDict
from A_modules.atomistic.file_processing.gromacs_file_parser.section_processing.section import (
    Section,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.base_handler import (
    BaseHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.default_handler import (
    DefaultHandler,
)
from typing import Dict, List, Optional, Tuple, OrderedDict
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.registry.handler_registry import (
    HandlerRegistry,
    handler_registry,
)


class FileSplitter:
    def __init__(self, handler_registry: HandlerRegistry = handler_registry):
        self.handler_registry = handler_registry
        self.suppressed_constructs: Optional[List[str]] = None

    # NOTE: make this more robust
    def split_file(self, filepath: str) -> OrderedDict[str, Section]:
        if filepath.endswith(".gro"):
            construct_type = "gro_file"
            handler_name = "gro"
        else:
            construct_type = None
            handler_name = "default"

        sections: OrderedDict[str, Section] = OrderedDict()
        current_section = Section(
            construct_type=construct_type, handler_name=handler_name
        )

        with open(filepath, "r") as file:
            for line in file:
                line = line.rstrip("\n")

                construct_type, name, handler_name = self._match_line(line)

                if construct_type:
                    key = self._generate_key(
                        sections, current_section.construct_type, current_section.name
                    )
                    sections[key] = current_section

                    current_section = Section(
                        construct_type=construct_type,
                        handler_name=handler_name,
                        name=name,
                    )

                current_section.add_line(line)

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
        base_key = f"{construct_type}_{name or 'no_name'}"

        if base_key in sections:
            last_key = next(reversed(sections), None)
            if last_key and last_key.startswith(base_key):
                return last_key
            raise ValueError(f"Duplicate section found: {base_key}")

        return base_key

    def _match_line(self, line: str) -> Tuple[str, str, str]:
        """Matches a line to a construct and returns its type, name, and handler name."""
        for handler_name, handler_class in self.handler_registry._handlers.items():
            if handler_class.re_pattern is None:
                continue

            match = handler_class.re_pattern.match(line)
            if match:
                name = match.group(1) if match.groups() else None
                self.suppressed_constructs = handler_class.suppress
                return handler_class.construct_name, name, handler_name

        # Fallback to DefaultHandler
        return None, None, DefaultHandler.construct_name
