import os
from typing import Callable, Dict


class FileCommenter:
    """
    Utility to add comments/notes to supported file types.
    """
    _comment_handlers: Dict[str, Callable] = {}

    def __init__(self):
        self._register_file_types()

    def _register_file_types(self):
        """
        Explicitly register file types and their handlers.
        """
        self._comment_handlers[".pdb"] = self._add_comment_pdb
        self._comment_handlers[".gro"] = self._add_comment_gro
        self._comment_handlers[".lmp"] = self._add_comment_lmp

    def add_comment(self, file_path: str, comment: str) -> None:
        """
        Add a comment to the top of a file using the appropriate handler.

        Args:
            file_path (str): Path to the file.
            comment (str): The comment to add.
        """
        file_extension = os.path.splitext(file_path)[1].lower()
        handler = self._comment_handlers.get(file_extension)

        if not handler:
            raise ValueError(f"Unsupported file type: {file_extension}. Supported types: {list(self._comment_handlers.keys())}")

        handler(file_path, comment)

    def _add_comment_pdb(self, file_path: str, comment: str) -> None:
        self._add_comment_to_file(file_path, comment, "#")

    def _add_comment_gro(self, file_path: str, comment: str) -> None:
        self._add_comment_to_file(file_path, comment, ";")

    def _add_comment_lmp(self, file_path: str, comment: str) -> None:
        self._add_comment_to_file(file_path, comment, "#")

    def _add_comment_to_file(self, file_path: str, comment: str, comment_symbol: str) -> None:
        """
        Add the comment to the top of the file.

        Args:
            file_path (str): The path to the file.
            comment (str): The comment to add.
            comment_symbol (str): The comment symbol for the file type.
        """
        with open(file_path, "r") as f:
            original_content = f.readlines()

        # Add the comment at the top
        formatted_comment = f"{comment_symbol} {comment.strip()}\n"
        updated_content = [formatted_comment] + original_content

        with open(file_path, "w") as f:
            f.writelines(updated_content)

    def supported_file_types(self):
        """
        Return all supported file types and their handlers.
        """
        return list(self._comment_handlers.keys())
