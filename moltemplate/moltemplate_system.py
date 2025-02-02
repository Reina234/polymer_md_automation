from moltemplate.polymer import MoltemplatePolymer
from moltemplate.solvent import MoltemplateSolvent
from modules.utils.shared.file_utils import check_directory_exists, check_file_type
from modules.rdkit.polymer_builders.base_polymer_generator import BasePolymerGenerator
from typing import List
import logging
import re

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class MoltemplateSystem:

    def __init__(
        self,
        n_units: int,
        polymer: BasePolymerGenerator,
        box_nm: List[float],
        openmscg_topol_path: str,
        cg_bond_length: float = 1.0,
        sol_resname: str = "SOL",
    ):
        check_file_type(openmscg_topol_path, "data")
        self.openmscg_topol_path = openmscg_topol_path
        self.open_mscg_mass_map = None
        self._process_openmscg_topol()
        self.box = box_nm
        self.polymer = MoltemplatePolymer(
            n_units=n_units,
            polymer_generator=polymer,
            cg_bond_length=cg_bond_length,
            open_mscg_mass_map=self.open_mscg_mass_map,
        )
        self.solvent = MoltemplateSolvent(
            open_mscg_mass_map=self.open_mscg_mass_map, sol_resname=sol_resname
        )

    def _process_openmscg_topol(self):
        mass_mapping = {}

        with open(self.openmscg_topol_path, "r") as file:
            lines = file.readlines()

        inside_masses = False
        for line in lines:
            line = line.strip()
            print(line)

            if line.strip() == "Masses":
                inside_masses = True
                continue

            if inside_masses:
                if line.strip().startswith("Atoms"):
                    break

                match = re.match(r"(\d+)\s+([\d.eE+-]+)\s+#\s*(\S+)", line)
                if match:
                    atom_type, mass, bead_name = match.groups()
                    mass_mapping[bead_name] = float(mass)

        self.open_mscg_mass_map = mass_mapping

    def _generate_masses_lt(self) -> str:
        mass_lines = ["Masses {"]

        for bead, mass in self.polymer.mass_mapping.items():
            mass_lines.append(f"  {bead} {mass}")

        mass_lines.append(f"  {self.solvent.resname} {self.solvent.mass}")
        mass_lines.append("}\n")
        return "\n".join(mass_lines)

    def write_system_lt(self, filename="system.lt", output_dir=None):
        if output_dir:
            check_directory_exists(output_dir, make_dirs=True)
            filename = f"{output_dir}/{filename}"

        with open(filename, "w") as f:
            f.write("# Moltemplate system file (with masses and full interactions)\n")
            f.write('import "forcefield.lt"\n\n')
            f.write(self._generate_masses_lt())
            f.write(self.polymer.generate_lt())
            f.write(self.solvent.generate_lt())

        logger.info(f"Moltemplate system file '{filename}' generated.")
