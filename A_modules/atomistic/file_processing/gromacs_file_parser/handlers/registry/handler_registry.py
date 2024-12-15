from typing import Dict, Type
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.base_handler import (
    BaseHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.base_handler import (
    BaseHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.section_handler import (
    SectionHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.includes_handler import (
    IncludesHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.conditional_if_handler import (
    ConditionalIfHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.default_handler import (
    DefaultHandler,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.gro_handler import (
    GroHandler,
)


class HandlerRegistry:
    def __init__(self):
        self._handlers: Dict[str, Type[BaseHandler]] = {}

    def register_handler(self, handler_class: Type[BaseHandler]) -> None:
        """Registers a handler class by its construct name."""
        construct_name = handler_class.construct_name
        if construct_name in self._handlers:
            raise ValueError(f"Handler '{construct_name}' is already registered.")
        self._handlers[construct_name] = handler_class

    def get_handler(self, handler_name: str) -> Type[BaseHandler]:
        """Retrieves a handler class by its name."""
        if handler_name not in self._handlers:
            available = ", ".join(self._handlers.keys())
            raise ValueError(
                f"Handler '{handler_name}' not found. Available handlers: {available}"
            )
        return self._handlers[handler_name]


# Initialize and populate the registry
handler_registry = HandlerRegistry()
handler_registry.register_handler(DefaultHandler)
handler_registry.register_handler(SectionHandler)
handler_registry.register_handler(IncludesHandler)
handler_registry.register_handler(ConditionalIfHandler)
handler_registry.register_handler(GroHandler)
