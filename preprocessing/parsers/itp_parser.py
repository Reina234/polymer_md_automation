import os
import logging
from typing import List, Dict

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class ITPParser:
    """
    A parser for GROMACS .itp files that handles extracting and modifying sections.
    """

    @staticmethod
    def read_file(file_path: str) -> List[str]:
        """
        Reads the content of a .itp file.

        Args:
            file_path (str): Path to the .itp file.

        Returns:
            List[str]: Lines from the file.
        """
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")

        with open(file_path, "r") as file:
            content = file.readlines()

        logger.info(f"[+] Successfully read file: {file_path}")
        return content

    @staticmethod
    def extract_section(content: List[str], section_name: str) -> List[str]:
        """
        Extracts a specific section from the .itp content.

        Args:
            content (List[str]): Lines from the .itp file.
            section_name (str): The name of the section to extract (e.g., "atomtypes").

        Returns:
            List[str]: Lines belonging to the specified section.
        """
        in_section = False
        section_lines = []

        for line in content:
            stripped = line.strip()
            if stripped.startswith(f"[ {section_name} ]"):
                in_section = True
                section_lines.append(line)  # Include the section header
                logger.info(f"[+] Found section: {section_name}")
                continue

            if in_section:
                if stripped.startswith("[") and not stripped.startswith(
                    f"[ {section_name} ]"
                ):
                    in_section = False
                    break
                section_lines.append(line)

        if not section_lines:
            logger.warning(f"[!] Section '{section_name}' not found in file.")
        return section_lines

    @staticmethod
    def remove_section(content: List[str], section_name: str) -> List[str]:
        """
        Removes a specific section from the .itp content.

        Args:
            content (List[str]): Lines from the .itp file.
            section_name (str): The name of the section to remove (e.g., "atomtypes").

        Returns:
            List[str]: Modified content with the section removed.
        """
        in_section = False
        updated_content = []

        for line in content:
            stripped = line.strip()
            if stripped.startswith(f"[ {section_name} ]"):
                in_section = True
                logger.info(f"[+] Removing section: {section_name}")
                continue

            if in_section:
                if stripped.startswith("[") and not stripped.startswith(
                    f"[ {section_name} ]"
                ):
                    in_section = False

            if not in_section:
                updated_content.append(line)

        return updated_content

    @staticmethod
    def append_section(
        base_content: List[str], section_content: List[str], section_name: str
    ) -> List[str]:
        """
        Appends a specific section to the end of a base .itp content.

        Args:
            base_content (List[str]): Lines of the base .itp file.
            section_content (List[str]): Lines of the section to append.
            section_name (str): Name of the section being appended (for logging).

        Returns:
            List[str]: Updated content with the section appended.
        """
        if not section_content:
            logger.warning(f"[!] Nothing to append for section '{section_name}'.")
            return base_content

        base_content.append("\n")  # Ensure a new line before appending
        base_content.extend(section_content)
        logger.info(f"[+] Appended section '{section_name}' to the base content.")
        return base_content

    @staticmethod
    def save_file(file_path: str, content: List[str]) -> None:
        """
        Saves the modified content to a .itp file.

        Args:
            file_path (str): Path to save the file.
            content (List[str]): Lines to write to the file.
        """
        with open(file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] Saved file: {file_path}")
