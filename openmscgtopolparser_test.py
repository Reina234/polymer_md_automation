from modules.lammps.parsers.open_mscg_data_parser import OpenMSCGDataParser

parser = OpenMSCGDataParser(data_file="cg_poly.data")

print(parser.bond_lengths)
