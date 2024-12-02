import os
import logging
from typing import List, Optional

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class TOPParser:
    """
    A parser specifically designed to read, restructure, and validate GROMACS .top files.
    """

    @staticmethod
    def read_file(input_file_path: str) -> List[str]:
        """
        Read the .top file and return its content as a list of lines.

        Args:
            input_file_path (str): Path to the .top file.

        Returns:
            List[str]: Lines from the .top file.
        """
        if not os.path.exists(input_file_path):
            raise FileNotFoundError(f"File not found: {input_file_path}")
        with open(input_file_path, "r") as file:
            content = file.readlines()
        logger.info(f"[+] .top file read successfully: {input_file_path}")
        return content

    @staticmethod
    def ensure_include_order(
        content: List[str],
        forcefield_include: str,
        monomer_include: str,
        solvent_include: Optional[str] = None,
    ) -> List[str]:
        """
        Ensure that includes are in the correct order:
        - Force field first.
        - Monomer parameters next.
        - Position restraints (if present) before solvent (if included).

        Args:
            content (List[str]): Lines from the .top file.
            forcefield_include (str): Path to the forcefield include.
            monomer_include (str): Path to the monomer topology include.
            solvent_include (Optional[str]): Path to the solvent topology include.

        Returns:
            List[str]: The modified content with includes in the correct order.
        """
        includes = {
            "forcefield": f'#include "{forcefield_include}"\n',
            "monomer": f'#include "{monomer_include}"\n',
        }
        if solvent_include:
            includes["solvent"] = f'#include "{solvent_include}"\n'

        posres_block = []
        in_posres = False

        # Filter out existing include lines and handle POSRES separately
        updated_content = []
        for line in content:
            stripped = line.strip()
            if stripped.startswith("#include"):
                # Skip existing includes
                continue
            elif stripped.startswith("#ifdef POSRES"):
                in_posres = True
                posres_block.append(line)
            elif in_posres:
                posres_block.append(line)
                if stripped.startswith("#endif"):
                    in_posres = False
            else:
                updated_content.append(line)

        # Insert includes in the correct order
        ordered_includes = [
            includes["forcefield"],
            includes["monomer"],
        ]
        if posres_block:
            ordered_includes += posres_block
        if "solvent" in includes:
            ordered_includes.append(includes["solvent"])

        updated_content = ordered_includes + updated_content

        logger.info("[+] Ensured correct order of includes.")
        return updated_content

    @staticmethod
    def handle_posres(content: List[str], posres: bool) -> List[str]:
        """
        Conditionally remove position restraints if POSRES is false.

        Args:
            content (List[str]): Lines from the .top file.
            posres (bool): Whether to retain position restraints.

        Returns:
            List[str]: The modified content with position restraints handled.
        """
        if posres:
            logger.info("[+] Retaining position restraints.")
            return content

        # Remove POSRES block
        in_posres = False
        updated_content = []
        for line in content:
            stripped = line.strip()
            if stripped.startswith("#ifdef POSRES"):
                in_posres = True
                logger.info("[+] Removing position restraints section.")
                continue
            if in_posres and stripped.startswith("#endif"):
                in_posres = False
                continue
            if not in_posres:
                updated_content.append(line)

        return updated_content

    @staticmethod
    def remove_section(content: List[str], section_name: str) -> List[str]:
        """
        Remove a specific section from the .top file.

        Args:
            content (List[str]): Lines from the .top file.
            section_name (str): Name of the section to remove (e.g., "defaults").

        Returns:
            List[str]: The modified content with the section removed.
        """
        in_section = False
        updated_content = []

        for line in content:
            stripped = line.strip()
            if stripped.startswith(f"[ {section_name} ]"):
                in_section = True
                logger.info(f"[+] Removing section: {section_name}")
                continue

            if in_section and stripped.startswith("[") and stripped.endswith("]"):
                in_section = False

            if not in_section:
                updated_content.append(line)

        return updated_content

    @staticmethod
    def save(output_file_path: str, content: List[str]) -> None:
        """
        Save the modified content back to a .top file.

        Args:
            output_file_path (str): Path to save the modified file.
            content (List[str]): The modified file content.
        """
        os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
        with open(output_file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] .top file saved to {output_file_path}")
