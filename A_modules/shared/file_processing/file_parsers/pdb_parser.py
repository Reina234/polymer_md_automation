import os
import logging
from typing import List, Dict, Optional
from pyparsing import (
    Word,
    nums,
    alphas,
    Combine,
    Literal,
    ParserElement,
    Optional as PPOptional,
    ParseException,
)

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


class PDBParser:
    """
    A parser for PDB files with modular grammar and parsing functionality.
    """

    class Grammar:
        """
        Encapsulates grammar definitions for parsing PDB files.
        """

        @staticmethod
        def atom_line() -> ParserElement:
            """
            Grammar for ATOM and HETATM lines.
            """
            record_name = Literal("ATOM") | Literal("HETATM")
            atom_name = Word(alphas + nums, max=4)
            residue_name = Word(alphas, max=3)
            x_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))
            y_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))
            z_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))

            # Define ATOM/HETATM line parsing
            return (
                record_name("record_name")
                + PPOptional(Word(nums))  # Atom serial (optional for parsing)
                + atom_name("atom_name")
                + PPOptional(Word(alphas, exact=1))  # Alternate location
                + residue_name("residue_name")
                + x_coord("x")
                + y_coord("y")
                + z_coord("z")
            )

        @staticmethod
        def cryst1_line() -> ParserElement:
            """
            Grammar for CRYST1 lines to extract box dimensions.
            """
            keyword = Literal("CRYST1")
            dimension = Combine(Word(nums) + PPOptional("." + Word(nums)))
            return keyword + dimension("x") + dimension("y") + dimension("z")

    def __init__(self):
        self.atom_line_grammar = self.Grammar.atom_line()
        self.cryst1_grammar = self.Grammar.cryst1_line()

    @staticmethod
    def read_file(input_file_path: str) -> List[str]:
        """
        Read the PDB file and return its content as a list of lines.

        Args:
            input_file_path (str): Path to the PDB file.

        Returns:
            List[str]: Lines from the PDB file.
        """
        if not os.path.exists(input_file_path):
            raise FileNotFoundError(f"File not found: {input_file_path}")
        with open(input_file_path, "r") as file:
            return file.readlines()

    def extract_atoms(self, content: List[str]) -> List[Dict[str, float]]:
        """
        Extract atom data from the PDB file content.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            List[Dict[str, float]]: List of atoms with their coordinates.
        """
        atoms = []
        for line in content:
            try:
                parsed = self.atom_line_grammar.parseString(line)
                atoms.append(
                    {
                        "atom_name": parsed.atom_name,
                        "residue_name": parsed.residue_name,
                        "x": float(parsed.x),
                        "y": float(parsed.y),
                        "z": float(parsed.z),
                    }
                )
            except ParseException:
                continue  # Skip lines that don't match
        logger.info(f"[+] Extracted {len(atoms)} atoms from PDB file.")
        return atoms

    def extract_box_dimensions(self, content: List[str]) -> Optional[List[float]]:
        """
        Extract box dimensions from the CRYST1 line in the PDB file.

        Args:
            content (List[str]): Lines from the PDB file.

        Returns:
            Optional[List[float]]: Box dimensions [x, y, z] in nm, or None if not found.
        """
        for line in content:
            try:
                parsed = self.cryst1_grammar.parseString(line)
                box_dimensions = [float(parsed.x), float(parsed.y), float(parsed.z)]
                logger.info(f"[+] Extracted box dimensions: {box_dimensions}")
                return box_dimensions
            except ParseException:
                continue
        logger.warning("[!] CRYST1 line not found; box dimensions are unavailable.")
        return None

    @staticmethod
    def save(output_file_path: str, content: List[str]) -> None:
        """
        Save the PDB content to a file.

        Args:
            output_file_path (str): Path to save the file.
            content (List[str]): PDB content to save.
        """
        os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
        with open(output_file_path, "w") as file:
            file.writelines(content)
        logger.info(f"[+] Saved PDB file to {output_file_path}")
