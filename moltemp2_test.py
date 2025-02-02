import numpy as np
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def extract_mass_mapping(file_path):
    mass_mapping = {}
    with open(file_path, "r") as file:
        lines = file.readlines()

    # Locate the "Masses" section
    inside_masses = False
    for line in lines:
        line = line.strip()
        if line.startswith("Masses"):
            inside_masses = True
            continue
        if inside_masses:
            if not line or line.startswith("Atoms"):  # Stop at the next section
                break
            match = re.match(r"(\d+)\s+([\d.]+)\s+#\s+(\S+)", line)
            if match:
                atom_type, mass, bead_name = match.groups()
                mass_mapping[bead_name] = float(mass)

    return mass_mapping


class Solvent:
    """Defines solvent molecules separately, with correct mass and parameters."""

    def __init__(self, sol_resname="SOL"):
        self.sol_resname = sol_resname
        self.mass = 86.0  # Solvent mass

    def generate_solvent_lt(self) -> str:
        """Creates Moltemplate definition for solvent."""
        return f"Molecule Solvent {{\n  $atom:{self.sol_resname} 0.0 0.0 0.0\n}}\n"


class Solvent:
    def __init__(self, sol_resname="SOL"):
        self.sol_resname = sol_resname

    def generate_solvent_lt(self):
        return (
            "Molecule Solvent {\n" f"  $atom:{self.sol_resname} 0.0 0.0 0.0\n" + "}\n"
        )


class Polymer:
    def __init__(
        self,
        N: int,
        polymer_generator: BasePolymerGenerator,
        cg_bond_length: float = 1.0,
    ):
        self.N = N
        self.polymer_generator = polymer_generator
        self.terminal_left = None
        self.terminal_right = None
        self.middle_beads = None
        self.cg_map = polymer_generator.cg_map
        self.n_repeats = None
        self.actual_n = None
        self.mass_mapping = {}  # Stores bead masses
        self._process_polymer_cg_map()
        self._generate_allowable_n()
        self._assign_masses()
        self.unique_bead_types = []
        self.cg_bond_length = cg_bond_length

    def process_polymer_cg_map(self) -> None:

        for bead in self.cg_map:
            bead_type = bead["bead_type"]
            mass = sum(bead["x-weight"])

            if bead_type not in self.unique_bead_types:
                self.unique_bead_types.append(bead_type)

            if bead_type not in self.mass_mapping:
                self.mass_mapping[bead_type] = mass
            else:
                stored_mass = self.mass_mapping[bead_type]

                if stored_mass == 0 and mass != 0:
                    self.mass_mapping[bead_type] = mass

                elif stored_mass != 0 and mass != 0 and abs(mass - stored_mass) > 1e-12:
                    logger.warning(
                        f"Bead type '{bead_type}' has conflicting masses: "
                        f"{stored_mass} vs {mass}. Retaining {stored_mass}"
                    )

        self.terminal_left = self.unique_bead_types[0]
        self.terminal_right = self.unique_bead_types[-1]
        self.middle_beads = self.unique_bead_types[1:-1]

    def _process_polymer_cg_map(self) -> None:
        """Extracts bead types and sets up sequence."""
        self.bead_mapping = self.polymer_generator.cg_map["bead_types"]
        bead_types = [bead["bead_type"] for bead in self.bead_mapping]
        bead_types = list(dict.fromkeys(bead_types))  # Remove duplicates
        self.terminal_left = bead_types[0]
        self.terminal_right = bead_types[-1]
        self.middle_beads = bead_types[1:-1]

    def _generate_allowable_n(self) -> None:
        """Determines the number of monomers allowed in the polymer."""
        repeating_units = len(self.middle_beads)
        self.n_repeats = (self.N - 2) // repeating_units
        self.actual_n = self.n_repeats * repeating_units + 2  # Actual polymer length

    def generate_polymer_lt(self) -> str:
        lines = ["Molecule Polymer {"]

        lines.append(f"  $atom:{self.terminal_left} 0.0 0.0 0.0")

        x_pos = 1.5
        for i in range(self.n_repeats):
            for bead in self.middle_beads:
                lines.append(f"  $atom:{bead}_{i+1} {x_pos:.1f} 0.0 0.0")
                x_pos += self.cg_bond_length

        lines.append(f"  $atom:{self.terminal_right} {x_pos:.1f} 0.0 0.0")

        # Add bonds
        lines.append("\n  # Bonds")
        all_atoms = (
            [self.terminal_left]
            + [f"{b}_{i+1}" for i in range(self.n_repeats) for b in self.middle_beads]
            + [self.terminal_right]
        )
        for i in range(len(all_atoms) - 1):
            lines.append(
                f"  bond:bond_{i+1} @atom:{all_atoms[i]} @atom:{all_atoms[i+1]}"
            )

        # Add angles
        lines.append("\n  # Angles")
        for i in range(len(all_atoms) - 2):
            lines.append(
                f"  angle:angle_{i+1} @atom:{all_atoms[i]} @atom:{all_atoms[i+1]} @atom:{all_atoms[i+2]}"
            )

        lines.append("}\n")
        return "\n".join(lines)

    def generate_mass_lt(self) -> str:
        mass_lines = ["Masses {"]
        for bead, mass in self.mass_mapping.items():
            mass_lines.append(f"  {bead} {mass}")
        mass_lines.append("}\n")
        return "\n".join(mass_lines)


class MoltemplateSystem:

    def __init__(
        self,
        N: int,
        polymer: BasePolymerGenerator,
        box,
        solvent_positions,
        ion_positions,
        sol_resname="SOL",
    ):
        self.N = N
        self.box = box
        self.solvent_positions = solvent_positions
        self.polymer = Polymer(N=N, polymer=polymer)
        self.solvent = Solvent(sol_resname=sol_resname)

    def write_system_lt(self, filename="system.lt"):
        """Writes the Moltemplate system file with polymer, solvent, and ions."""
        with open(filename, "w") as f:
            f.write("# Moltemplate system file (generated by Python)\n")
            f.write('import "forcefield.lt"\n\n')
            f.write(self.polymer.generate_polymer_lt())
            f.write(self.solvent.generate_solvent_lt())
            f.write(self.ion.generate_ion_lt())

            # Instantiate polymer in center of box
            x_center = 0.5 * (self.box[0] + self.box[1])
            y_center = 0.5 * (self.box[2] + self.box[3])
            z_center = 0.5 * (self.box[4] + self.box[5])
            f.write(
                f"new Polymer polymer_instance {{ translate {x_center} {y_center} {z_center} }}\n\n"
            )

            # Instantiate solvent
            for i, pos in enumerate(self.solvent_positions):
                f.write(
                    f"new Solvent solvent_{i+1} {{ translate {pos[0]} {pos[1]} {pos[2]} }}\n"
                )

            # Instantiate ions
            for j, pos in enumerate(self.ion_positions):
                f.write(
                    f"new Ion ion_{j+1} {{ translate {pos[0]} {pos[1]} {pos[2]} }}\n"
                )

        print(f"Moltemplate system file '{filename}' generated.")


# --------------------------
# Support Functions for Bead Placement
# --------------------------
def create_grid(box, spacing):
    """Creates a regular grid inside the box."""
    x_min, x_max, y_min, y_max, z_min, z_max = box
    x_coords = np.arange(x_min, x_max, spacing)
    y_coords = np.arange(y_min, y_max, spacing)
    z_coords = np.arange(z_min, z_max, spacing)
    return [(x, y, z) for x in x_coords for y in y_coords for z in z_coords]


def remove_overlaps(solvent_positions, polymer_positions, cutoff):
    """Removes solvent beads too close to polymer beads."""
    filtered = [
        s
        for s in solvent_positions
        if all(
            np.linalg.norm(np.array(s) - np.array(p)) >= cutoff
            for p in polymer_positions
        )
    ]
    return filtered


# --------------------------
# Moltemplate System Builder
# --------------------------
class MoltemplateSystem:
    """Generates a complete Moltemplate system including polymer, solvent, and ions."""

    def __init__(
        self, N, box, solvent_positions, polymer_charge, ion_spacing, polymer_generator
    ):
        self.box = box
        self.solvent_positions = solvent_positions
        self.polymer = Polymer(N, polymer_generator)
        self.solvent = Solvent()
        self.ion = Ion(polymer_charge, box, ion_spacing)

    def write_system_lt(self, filename="system.lt"):
        """Writes the Moltemplate system file."""
        with open(filename, "w") as f:
            f.write("# Moltemplate system file (auto-neutralized)\n")
            f.write('import "forcefield.lt"\n\n')
            f.write(self.polymer.generate_polymer_lt())
            f.write(self.solvent.generate_solvent_lt())
            f.write(self.ion.generate_ion_lt())

            # Instantiate molecules
            x_center, y_center, z_center = [
                (self.box[i] + self.box[i + 1]) / 2 for i in range(0, 6, 2)
            ]
            f.write(
                f"new Polymer polymer_instance {{ translate {x_center} {y_center} {z_center} }}\n\n"
            )

            for i, pos in enumerate(self.solvent_positions):
                f.write(
                    f"new Solvent solvent_{i+1} {{ translate {pos[0]} {pos[1]} {pos[2]} }}\n"
                )

            for j, pos in enumerate(self.ion.ion_positions):
                f.write(
                    f"new Ion ion_{j+1} {{ translate {pos[0]} {pos[1]} {pos[2]} }}\n"
                )

        print(f"Moltemplate system file '{filename}' generated.")


# --------------------------
# LAMMPS Input File Generator
# --------------------------
class LAMMPSInputGenerator:
    """Generates input scripts for LAMMPS minimization and equilibration runs."""

    def generate_minimization_input(self, data_file):
        """Creates an energy minimization script."""
        with open("in.minim", "w") as f:
            f.write(
                f"""units real
atom_style full
read_data {data_file}

min_style cg
minimize 1.0e-4 1.0e-6 100 1000
"""
            )

    def generate_equilibration_inputs(self, data_file):
        """Creates NVT and NPT equilibration scripts."""
        with open("in.nvt", "w") as f:
            f.write(
                f"""units real
atom_style full
read_data {data_file}

fix 1 all nvt temp 300.0 300.0 100.0
run 50000
"""
            )

        with open("in.npt", "w") as f:
            f.write(
                f"""units real
atom_style full
read_data {data_file}

fix 1 all npt temp 300.0 300.0 100.0 iso 1.0 1.0 1000.0
run 50000
"""
            )


# --------------------------
# Main Driver Function
# --------------------------
def main():
    N = 10  # Number of polymer beads
    box = (-20.0, 20.0, -20.0, 20.0, -20.0, 20.0)
    solvent_spacing = 2.0
    ion_spacing = 5.0
    overlap_cutoff = 1.5
    polymer_positions = [(i, 0, 0) for i in range(N)]

    solvent_positions = remove_overlaps(
        create_grid(box, solvent_spacing), polymer_positions, overlap_cutoff
    )
    ion_positions = create_grid(box, ion_spacing)[:5]

    system = MoltemplateSystem(N, box, solvent_positions, ion_positions)
    system.write_system_lt("system.lt")


if __name__ == "__main__":
    main()
