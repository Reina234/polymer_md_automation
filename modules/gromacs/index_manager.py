import os
import subprocess
import logging
from typing import List
from modules.utils.shared.file_utils import check_directory_exists
from config.paths import TEMP_DIR

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GromacsIndexManager:
    index_file = "index.ndx"

    def __init__(
        self,
        gro_file: str,
        output_dir: str = TEMP_DIR,
        poly_name: str = "UNL",
        ion_names: List[str] = None,
    ):
        if ion_names is None:
            ion_names = ["NA", "CL"]
        self.gro_file = gro_file
        check_directory_exists(output_dir, make_dirs=True)
        self.index_file_path = os.path.join(output_dir, self.index_file)
        self.poly_name = poly_name
        self.ion_names = ion_names

        cmin, cmax = self.parse_gro_for_polymer_carbons(
            self.gro_file, polymer_resname=self.poly_name
        )
        self.first_carbon = cmin
        self.last_carbon = cmax

    def _filter_existing_ions(self) -> List[str]:

        with open(self.gro_file, "r") as f:
            gro_contents = f.read()

        existing = [ion for ion in self.ion_names if f" {ion} " in gro_contents]
        if not existing:
            logger.info("No ions found in the GRO file.")
        return existing

    def create_index_file(self) -> str:
        logger.info(f"Creating index file at: {self.index_file_path}")

        existing_ions = self._filter_existing_ions()

        instructions = []
        next_group_num = 4

        instructions.append(f"r {self.poly_name}")
        instructions.append(f"name {next_group_num} Polymer")
        polymer_group_num = next_group_num
        next_group_num += 1

        if existing_ions:
            ion_string = " ".join(existing_ions)
            instructions.append(f"r {ion_string}")
            instructions.append(f"name {next_group_num} Ions")
            ions_group_num = next_group_num
            next_group_num += 1

            instructions.append(f"0 & !{polymer_group_num} & !{ions_group_num}")
            instructions.append(f"name {next_group_num} Solvent")
            next_group_num += 1
        else:

            instructions.append(f"0 & !{polymer_group_num}")
            instructions.append(f"name {next_group_num} Solvent")
            next_group_num += 1

        if self.first_carbon is not None and self.last_carbon is not None:

            instructions.append(
                f"a {self.first_carbon} |a {self.last_carbon}  & r {self.poly_name}"
            )
            instructions.append(f"name {next_group_num} Polymer_Carbons_Start_and_End")
            next_group_num += 1

        else:
            logger.warning(
                "No polymer carbons found in .gro; skipping start/end carbon groups."
            )

        instructions.append("l")
        instructions.append("q")

        ndx_script = "\n".join(instructions) + "\n"

        cmd = ["gmx", "make_ndx", "-f", self.gro_file, "-o", self.index_file_path]
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate(input=ndx_script)
        p.wait()

        if not os.path.exists(self.index_file_path):
            raise FileNotFoundError(
                f"Failed to create index file: {self.index_file_path}"
            )

        logger.info(f"Index file created at {self.index_file_path}")
        return self.index_file_path

    @staticmethod
    def parse_gro_for_polymer_carbons(gro_file: str, polymer_resname: str = "UNL"):
        if not os.path.exists(gro_file):
            raise FileNotFoundError(f"GRO file not found: {gro_file}")

        with open(gro_file, "r") as f:
            lines = f.readlines()

        if len(lines) < 3:
            return None, None, None, None

        try:
            atom_count = int(lines[1].strip())
        except ValueError:
            atom_count = 0

        atom_lines = lines[2 : 2 + atom_count]

        min_carbon_index = None
        min_carbon_name = None
        max_carbon_index = None
        max_carbon_name = None

        for line in atom_lines:
            line_str = line.rstrip("\n")

            residue_name = line_str[5:10].strip()
            atom_name = line_str[10:15].strip()
            try:
                atom_index = int(line_str[15:20])
            except ValueError:
                continue

            if residue_name == polymer_resname and atom_name.upper().startswith("C"):
                if min_carbon_index is None or atom_index < min_carbon_index:
                    min_carbon_index = atom_index
                    min_carbon_name = atom_name

                if max_carbon_index is None or atom_index > max_carbon_index:
                    max_carbon_index = atom_index
                    max_carbon_name = atom_name

        return min_carbon_name, max_carbon_name
