from A_modules.atomistic.file_processing.gromacs_file_parser.file_splitter import (
    FileSplitter,
)
from A_modules.atomistic.file_processing.gromacs_file_parser.handlers.registry.handler_registry import (
    handler_registry,
)

file_splitter = FileSplitter()
sections = file_splitter.split_file("TEST/POLY_GMX.gro")
print(sections)

test_section = next(iter(sections.values()))
# print(test_section.name)
print(test_section.lines)

tester = handler_registry.get_handler(test_section.handler_name)()
tester.process(test_section)
print(tester.content)

test2 = tester.export()
print(test2.lines)
