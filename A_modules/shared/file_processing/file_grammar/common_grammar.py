from pyparsing import Word, alphas, nums, Combine, LineStart, Optional, Group


class CommonGrammar:
    integer = Word(nums).setParseAction(lambda t: int(t[0]))
    real = Combine(Word(nums) + "." + Word(nums)).setParseAction(lambda t: float(t[0]))
    identifier = Word(alphas)
    newline = LineStart()
