from A_modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser

parser = GromacsParser()
sections = parser.parse("TRIAL/acpype_output/POLY_GMX.top")
print(sections)
