from ZZZ_parser.file_splitter import FileSplitter
from ZZZ_parser.handler_registry import handler_registry

file_splitter = FileSplitter()
sections = file_splitter.split_file("hexane.top")
print(sections)

for key, section in sections.items():
    print(f"Section Key: {key}")
    print(f"Section Name: {section.name}")
    print(f"Section Lines: {section.lines}")  # Actual content of the lines
    print(f"Number of Lines: {len(section.lines)}")  # Count of lines


test_section = sections["section_molecules"]
print(test_section.name)
print(test_section.lines)

tester = handler_registry.get_handler(test_section.handler_name)()
tester.process(test_section)
print(tester.content)
