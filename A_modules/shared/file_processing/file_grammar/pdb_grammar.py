from pyparsing import (
    Literal,
    Word,
    nums,
    alphas,
    Combine,
    Suppress,
    LineEnd,
    Optional as PPOptional,
    ParseResults,
    ParserElement,
)


class PDBGrammar:
    """
    Defines pyparsing grammars for parsing PDB files.
    """

    def __init__(self):
        self.atom_line_grammar = self._define_atom_line_grammar()
        self.cryst1_grammar = self._define_cryst1_grammar()

    @staticmethod
    def _define_atom_line_grammar() -> ParserElement:
        """
        Define the grammar for parsing ATOM and HETATM lines.
        """
        record_name = Literal("ATOM") | Literal("HETATM")
        atom_serial = Word(nums)
        atom_name = Word(alphas + nums, exact=4).setParseAction(lambda t: t[0].strip())
        residue_name = Word(alphas, max=3)
        x_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))
        y_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))
        z_coord = Combine(Word(nums) + PPOptional("." + Word(nums)))

        # Fixed-width grammar for ATOM and HETATM lines
        atom_line = (
            record_name("record_name")
            + Suppress(Word(nums))  # Skip serial number (irrelevant for this example)
            + atom_name("atom_name")
            + Suppress(Word(alphas, exact=1))  # Skip alternative location indicator
            + residue_name("residue_name")
            + Suppress(Word(nums))  # Skip residue sequence number
            + x_coord("x")
            + y_coord("y")
            + z_coord("z")
        )
        return atom_line

    @staticmethod
    def _define_cryst1_grammar() -> ParserElement:
        """
        Define the grammar for parsing CRYST1 lines.
        """
        keyword = Literal("CRYST1")
        dimension = Combine(Word(nums) + PPOptional("." + Word(nums)))
        cryst1_line = (
            keyword
            + dimension("x")
            + dimension("y")
            + dimension("z")
            + Suppress(Word(nums))  # Skip angles
        )
        return cryst1_line
