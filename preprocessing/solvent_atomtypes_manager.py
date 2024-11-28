import os
import logging
from typing import List

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class SolventAtomtypesManager:
    """
    Manages solvent atomtypes, including extraction, storage, and retrieval.
    """

    ATOMTYPES_DIR = "solvent_atomtypes"

    def __init__(self):
        os.makedirs(self.ATOMTYPES_DIR, exist_ok=True)

    def extract_and_store_atomtypes(self, solvent_file: str, solvent_name: str) -> None:
        """
        Extracts the [ atomtypes ] section from a solvent file and stores it.

        Args:
            solvent_file (str): Path to the solvent .itp file.
            solvent_name (str): Name of the solvent to use as a key.
        """
        content = self._read_file(solvent_file)
        atomtypes_section = self._extract_section(content, "atomtypes")
        if not atomtypes_section:
            raise ValueError(f"No [ atomtypes ] section found in {solvent_file}.")

        # Save the atomtypes to the solvent_atomtypes directory
        atomtypes_path = os.path.join(self.ATOMTYPES_DIR, f"{solvent_name}.itp")
        self._save_file(atomtypes_path, atomtypes_section)
        logger.info(f"[+] Saved atomtypes for {solvent_name} to {atomtypes_path}.")

        # Remove the atomtypes section from the original file
        updated_content = self._remove_section(content, "atomtypes")
        self._save_file(solvent_file, updated_content)
        logger.info(f"[+] Removed [ atomtypes ] section from {solvent_file}.")

    def retrieve_atomtypes(self, solvent_name: str) -> List[str]:
        """
        Retrieves the [ atomtypes ] section for a specified solvent.

        Args:
            solvent_name (str): Name of the solvent.

        Returns:
            List[str]: Lines of the [ atomtypes ] section.
        """
        atomtypes_path = os.path.join(self.ATOMTYPES_DIR, f"{solvent_name}.itp")
        if not os.path.exists(atomtypes_path):
            raise FileNotFoundError(f"Atomtypes for {solvent_name} not found.")
        return self._read_file(atomtypes_path)

    def _read_file(self, file_path: str) -> List[str]:
        """Helper to read a file."""
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"File not found: {file_path}")
        with open(file_path, "r") as file:
            return file.readlines()

    def _save_file(self, file_path: str, content: List[str]) -> None:
        """Helper to save content to a file."""
        with open(file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] Saved file: {file_path}")

    def _extract_section(self, content: List[str], section_name: str) -> List[str]:
        """Helper to extract a specific section."""
        in_section = False
        section_lines = []
        for line in content:
            if line.strip().startswith(f"[ {section_name} ]"):
                in_section = True
                section_lines.append(line)
                continue
            if in_section:
                if line.strip().startswith("[") and not line.strip().startswith(
                    f"[ {section_name} ]"
                ):
                    break
                section_lines.append(line)
        return section_lines

    def _remove_section(self, content: List[str], section_name: str) -> List[str]:
        """Helper to remove a specific section."""
        in_section = False
        updated_content = []
        for line in content:
            if line.strip().startswith(f"[ {section_name} ]"):
                in_section = True
                logger.info(f"[+] Removing section: {section_name}")
                continue
            if in_section:
                if line.strip().startswith("[") and not line.strip().startswith(
                    f"[ {section_name} ]"
                ):
                    in_section = False
            if not in_section:
                updated_content.append(line)
        return updated_content
