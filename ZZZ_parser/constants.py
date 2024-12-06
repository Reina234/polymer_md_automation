import re
from typing import Dict, List, Optional, Type, Union
from dataclasses import dataclass
from ZZZ_parser.handlers.base_handler import BaseHandler
from ZZZ_parser.handlers.section_handler import SectionHandler
from ZZZ_parser.handlers.includes_handler import IncludesHandler
from ZZZ_parser.handlers.conditional_if_handler import ConditionalIfHandler


CONDITIONAL_HANDLER = "ConditionalHandler"
INCLUDE_HANDLER = "IncludeHandler"
SECTION_HANDLER = "SectionHandler"
DEFAULT_HANDLER = "DefaultHandler"


HANDLERS: Dict[str, BaseHandler] = {
    CONDITIONAL_HANDLER: ConditionalIfHandler,
    INCLUDE_HANDLER: IncludesHandler,
    SECTION_HANDLER: SectionHandler,
}


def get_handler(
    handler_name: str, handler_list: Dict[str, BaseHandler] = HANDLERS
) -> BaseHandler:
    if handler_name in handler_list:
        return handler_list[handler_name]()
    else:
        raise ValueError(
            f"Handler '{handler_name}' not found in handler list. Available handlers: {list(handler_list.keys())}"
        )


@dataclass
class GROMACSConstructs:
    pattern: re.Pattern
    handler: str
    suppress: Optional[List[str]]


GROMACS_CONSTRUCTS = {
    "section": GROMACSConstructs(
        pattern=re.compile(r"^\s*\[\s*(.+?)\s*\]\s*$"),
        handler=SECTION_HANDLER,
        suppress=None,
    ),
    "include": GROMACSConstructs(
        pattern=re.compile(r'^\s*#\s*include\s+"(.+?)"\s*$'),
        handler=INCLUDE_HANDLER,
        suppress=None,
    ),
    "conditional_if": GROMACSConstructs(
        pattern=re.compile(r"^\s*#\s*(ifdef|ifndef)\s+.*$"),
        handler=CONDITIONAL_HANDLER,
        suppress=["include"],
    ),
}
