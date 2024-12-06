from ZZZ_parser.file_splitter import FileSplitter

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

from ZZZ_parser.constants import HANDLERS

handler = HANDLERS[test_section.handler_type]()
handler.process(test_section)

print(handler.content)
print(handler.top_line)
tester2_section = handler.export()

from ZZZ_parser.constants import get_handler_new

# NOTE: get handler might be more readable
handler = get_handler_new(test_section.handler_type)
handler.process(test_section)

print(handler.content)
print(handler.top_line)
tester2_section = handler.export()
