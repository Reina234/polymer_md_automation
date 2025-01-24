from mscg.cli import cgdump

cgdump.main(file="result.p", dump=["Pair_B2-B2,0.05,10.0,0.05,L2"])
import matplotlib

matplotlib.use("Agg")  # Use non-GUI backend

import matplotlib.pyplot as plt
import numpy as np
from mscg.cli import cgdump

# Extract table from MSCG data
cgdump.main(file="result.p", dump=["Pair_B2-B2,0.05,10.0,0.05,L2"])

# Load table
table = np.loadtxt("Pair_B2-B2.table", skiprows=5)

# Plot forces and potential energy
plt.plot(table[:, 1], table[:, 3], label="Force - Kcal/mol/angstrom")
plt.plot(table[:, 1], table[:, 2], label="Potential Energy - Kcal/mol")

plt.legend(loc="upper right")
plt.xlim(2, 10)
plt.ylim(-2, 6)

# Save the plot instead of showing it
plt.savefig("force_potential_plot.png", dpi=300)
print("Plot saved as force_potential_plot.png. Open the file to view the plot.")
