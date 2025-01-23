import MDAnalysis as mda

# Load trajectory
u = mda.Universe("temp/npt_short_1.tpr", "temp/npt_short_1.trr")

# Loop through frames and check forces
for ts in u.trajectory[:5]:  # Check first 5 frames
    forces = u.atoms.forces
    print(f"Frame {ts.frame}: First 5 Forces: {forces[:5]}")

    # If all zeros or NaN, FM will fail!
    if (forces == 0).all():
        print(f"[ERROR] Frame {ts.frame}: All forces are zero!")
