from A_modules.atomistic.gromacs.parser.gromacs_parser import (
    GromacsParser,
)
from A_modules.atomistic.gromacs.parser.registries.handler_registry import (
    handler_registry,
)

file_splitter = GromacsParser()
sections = file_splitter.parse("TEST/POLY_GMX.gro")
print(sections)

test_section = next(iter(sections.values()))
# print(test_section.name)
print(test_section.lines)

tester = handler_registry.get_handler(test_section.handler_name)()
tester.process(test_section)
print(tester.content)

test2 = tester.export()
print(test2.lines)
print(tester.box_dimensions)
