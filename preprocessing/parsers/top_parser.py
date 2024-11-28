import os
import logging
from typing import List, Optional

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


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
            # Detect section start
            if stripped.startswith(f"[ {section_name} ]"):
                in_section = True
                logger.info(f"[+] Removing section: {section_name}")
                continue

            # Exit the section when encountering a new header
            if in_section and stripped.startswith("[") and stripped.endswith("]"):
                in_section = False

            if not in_section:
                updated_content.append(line)

        return updated_content

    @staticmethod
    def ensure_include_order(
        content: List[str], includes: List[str], solvent_include: str
    ) -> List[str]:
        """
        Ensure that includes are in the correct order:
        - Force fields before .itp files
        - Solvent topology after position restraints, if any

        Args:
            content (List[str]): Lines from the .top file.
            includes (List[str]): List of include paths to ensure.
            solvent_include (str): The include path for the solvent topology.

        Returns:
            List[str]: The modified content with includes in the correct order.
        """
        # Separate includes into force fields, .itp files, and solvent
        force_field_includes = [inc for inc in includes if inc.endswith(".ff")]
        itp_includes = [inc for inc in includes if inc.endswith(".itp")]
        solvent_line = f'#include "{solvent_include}"\n'

        # Remove all includes from the content
        content = [line for line in content if not line.strip().startswith("#include")]

        # Add force field includes first
        for include in force_field_includes:
            content.insert(0, f'#include "{include}"\n')

        # Add .itp files next
        for include in itp_includes:
            content.append(f'#include "{include}"\n')

        # Handle solvent include placement
        for i, line in enumerate(content):
            if line.strip().startswith("#ifdef POSRES"):
                # Add the solvent include after the position restraint block
                for j in range(i + 1, len(content)):
                    if content[j].strip().startswith("#endif"):
                        content.insert(j + 1, solvent_line)
                        break
                break
        else:
            # Add solvent include at the end if no POSRES block
            content.append(solvent_line)

        logger.info("[+] Ensured correct order of includes.")
        return content

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
    def validate_structure(content: List[str], required_sections: List[str]) -> bool:
        """
        Validate that the required sections are present in the correct order.

        Args:
            content (List[str]): Lines from the .top file.
            required_sections (List[str]): List of required section names.

        Returns:
            bool: True if the structure is valid, False otherwise.
        """
        section_order = []
        for line in content:
            stripped = line.strip()
            if stripped.startswith("[") and stripped.endswith("]"):
                section_name = stripped[1:-1].strip()
                section_order.append(section_name)

        # Validate the order of sections
        for section in required_sections:
            if section not in section_order:
                logger.error(f"[!] Missing required section: {section}")
                return False
        logger.info("[+] File structure is valid.")
        return True

    @staticmethod
    def save(output_file_path: str, content: List[str]) -> None:
        """
        Save the modified content back to a .top file.

        Args:
            output_file_path (str): Path to save the modified file.
            content (List[str]): The modified file content.
        """
        with open(output_file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] .top file saved to {output_file_path}")
