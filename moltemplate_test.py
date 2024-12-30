import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array


import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


def add_polymer_to_solvent(polymer_file, solvent_file, output_file, cutoff=0.2):
    """
    Add a polymer to a pre-equilibrated solvent box, removing overlapping solvent molecules.

    :param polymer_file: Path to the polymer GRO file.
    :param solvent_file: Path to the solvent GRO file.
    :param output_file: Path to save the combined system.
    :param cutoff: Distance cutoff (in nm) for removing overlapping solvent molecules.
    """
    u_polymer = mda.Universe(polymer_file)
    u_solvent = mda.Universe(solvent_file)

    # Center polymer in the solvent box
    polymer_center = u_polymer.atoms.center_of_mass()
    solvent_center = u_solvent.dimensions[:3] / 2
    u_polymer.atoms.translate(solvent_center - polymer_center)

    # Identify overlapping solvent molecules
    distances = distance_array(
        u_polymer.atoms.positions, u_solvent.atoms.positions, box=u_solvent.dimensions
    )
    overlapping_solvent_indices = distances.min(axis=0) < cutoff

    # Create a selection of non-overlapping solvent molecules
    solvent_atoms = u_solvent.atoms
    non_overlapping_solvent = solvent_atoms[~overlapping_solvent_indices]

    # Combine polymer and non-overlapping solvent
    combined = mda.Merge(u_polymer.atoms, non_overlapping_solvent)
    combined.dimensions = u_solvent.dimensions  # Retain original box dimensions

    # Save the combined system
    with mda.Writer(output_file, n_atoms=combined.atoms.n_atoms) as W:
        W.write(combined)


# Usage
add_polymer_to_solvent(
    "trial.gro",
    "TEST_RUN_28/hexane_run_1/equilibriated_gros/temp_298.gro",
    "output.gro",
    cutoff=0.2,
)
