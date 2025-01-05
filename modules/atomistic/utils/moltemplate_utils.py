import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from modules.shared.utils.file_utils import prepare_output_file_path
from typing import Optional


def add_polymer_to_solvent(
    polymer_file: str,
    solvent_file: str,
    output_dir: Optional[str] = None,
    output_name: Optional[str] = None,
    cutoff=0.2,
) -> str:
    """
    Add a polymer to a pre-equilibrated solvent box, removing overlapping solvent molecules.

    :param polymer_file: Path to the polymer GRO file.
    :param solvent_file: Path to the solvent GRO file.
    :param output_dir: Directory to save the combined system.
    :param output_name: Name of the output file.
    :param cutoff: Distance cutoff (in nm) for removing overlapping solvent molecules.
    """
    output_path = prepare_output_file_path(polymer_file, "gro", output_dir, output_name)
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
    overlapping_atoms = distances.min(axis=0) < cutoff
    overlapping_residues = u_solvent.atoms[overlapping_atoms].residues.resids

    # Select non-overlapping solvent molecules
    non_overlapping_solvent = u_solvent.select_atoms(
        f"not resid {' '.join(map(str, overlapping_residues))}"
    )

    # Combine polymer and non-overlapping solvent
    combined = mda.Merge(u_polymer.atoms, non_overlapping_solvent)
    combined.dimensions = u_solvent.dimensions  # Retain original box dimensions

    # Save the combined system
    with mda.Writer(output_path, n_atoms=combined.atoms.n_atoms) as W:
        W.write(combined)

    return output_path
