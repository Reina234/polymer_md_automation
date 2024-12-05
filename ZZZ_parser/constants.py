import re
from typing import Dict, List, Optional, Type
from dataclasses import dataclass
from ZZZ_parser.handlers.base_handler import BaseHandler
from ZZZ_parser.handlers.section_handler import SectionHandler
from ZZZ_parser.handlers.includes_handler import IncludesHandler
from ZZZ_parser.handlers.conditional_if_handler import ConditionalIfHandler


CONDITIONAL_HANDLER = "ConditionalHandler"
INCLUDE_HANDLER = "IncludeHandler"
SECTION_HANDLER = "SectionHandler"
DEFAULT_HANDLER = "DefaultHandler"


HANDLERS = {
    CONDITIONAL_HANDLER: ConditionalIfHandler,
    INCLUDE_HANDLER: IncludesHandler,
    SECTION_HANDLER: SectionHandler,
}


@dataclass
class GROMACSConstructs:
    pattern: re.Pattern
    handler: Type[BaseHandler]
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
