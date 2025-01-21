import MDAnalysis as mda
import networkx as nx
from collections import defaultdict

# from hoomd import data
# import gsd.hoomd

# Load the GROMACS atomistic trajectory
u = mda.Universe(
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.tpr",
    "test_29_full/styrene/gromacs_output/temp_298_compressibility_0.000124.xtc",
)


# Create graph to detect repeat units
G = nx.Graph()
for bond in u.bonds:
    G.add_edge(bond[0], bond[1])

# Detect repeat units (connected components)
repeat_units = list(nx.connected_components(G))

# Map each repeat unit to a CG bead
mapping = defaultdict(list)
for i, unit in enumerate(repeat_units):
    bead_name = f"B{i+1}"  # B1, B2, ...
    mapping[bead_name] = [atom.name for atom in unit]

# Write mapping file
with open("mapping.txt", "w") as f:
    f.write("[PVP]\n")
    for bead, atoms in mapping.items():
        f.write(f"{bead} P3 {' '.join(atoms)}\n")
