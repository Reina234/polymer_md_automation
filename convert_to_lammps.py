import MDAnalysis as mda

# Load the XTC trajectory with its corresponding GRO file
u = mda.Universe("out.gro", "out.xtc")


# Extract box dimensions from the GRO file
box_dims = u.dimensions  # [x_min, x_max, y_min, y_max, z_min, z_max]

# Open LAMMPS trajectory file for writing
with open("cg.lammpstrj", "w") as lmp:
    for ts in u.trajectory:
        num_atoms = len(u.atoms)

        # Ensure all box dimensions are valid
        x_min, x_max = 0.0, box_dims[0] if box_dims[0] > 0 else 10.0
        y_min, y_max = 0.0, box_dims[1] if box_dims[1] > 0 else 10.0
        z_min, z_max = 0.0, box_dims[2] if box_dims[2] > 0 else 10.0

        # Write LAMMPS headers
        lmp.write(f"ITEM: TIMESTEP\n{ts.frame}\n")
        lmp.write(f"ITEM: NUMBER OF ATOMS\n{num_atoms}\n")
        lmp.write("ITEM: BOX BOUNDS pp pp pp\n")
        lmp.write(f"{x_min} {x_max}\n")
        lmp.write(f"{y_min} {y_max}\n")
        lmp.write(f"{z_min} {z_max}\n")
        lmp.write("ITEM: ATOMS id type x y z\n")

        # Write atom positions
        for i, atom in enumerate(u.atoms):
            lmp.write(
                f"{i+1} {atom.type} {atom.position[0]} {atom.position[1]} {atom.position[2]}\n"
            )

print("Conversion to cg.lammpstrj completed successfully!")
