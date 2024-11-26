import os
from typing import Callable, Dict


class FileCommenter:
    """
    Utility to add comments/notes to supported file types with enforced registration via decorators.
    """
    _comment_handlers: Dict[str, Callable] = {}

    @classmethod
    def register_file_type(cls, extension: str):
        """
        Decorator to register a file type and its associated comment handler.

        Args:
            extension (str): File extension (e.g., ".pdb").
        """
        def decorator(func: Callable):
            if not extension.startswith("."):
                raise ValueError(f"File extension must start with a dot (e.g., '.pdb'). Got: {extension}")
            cls._comment_handlers[extension.lower()] = func
            return func
        return decorator

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

        handler(self, file_path, comment)

    @register_file_type(".pdb")
    def _add_comment_pdb(self, file_path: str, comment: str) -> None:
        self._add_comment_to_file(file_path, comment, "#")

    @register_file_type(".gro")
    def _add_comment_gro(self, file_path: str, comment: str) -> None:
        self._add_comment_to_file(file_path, comment, ";")

    @register_file_type(".lmp")
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

    @classmethod
    def supported_file_types(cls):
        """
        Return all supported file types and their handlers.
        """
        return list(cls._comment_handlers.keys())


# Example Usage
if __name__ == "__main__":
    commenter = FileCommenter()

    # Add a comment to a PDB file
    commenter.add_comment("example.pdb", "This is a PDB file.")

    # List supported file types
    print("Supported file types:", FileCommenter.supported_file_types())
