from typing import Optional


handlers = {}


class Section:
    def __init__(
        self, construct_type: str, handler_type: str, name: Optional[str] = None
    ):
        self.construct_type = construct_type
        self.handler_type = handler_type  # The handler type, e.g., "SectionHandler"
        self.name = name  # The name of the section, extracted from the match
        self.lines = []  # Lines belonging to this section

    def add_line(self, line: str):
        """
        Adds a line to the section's content.

        :param line: The line to add.
        """
        self.lines.append(line)

    def __repr__(self):
        """
        String representation of the Section for debugging purposes.
        """
        return f"Section(construct_type={self.construct_type}, handler={self.handler_type}, name={self.name}, lines={len(self.lines)})"
