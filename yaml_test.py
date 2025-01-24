import yaml
import pandas as pd
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.gro_handler import GroHandler
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from typing import Optional, Dict, List
from modules.atomistic.gromacs.parser.itp_parser import ITPParser

import yaml
import pandas as pd
from modules.atomistic.gromacs.parser.itp_parser import ITPParser
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from typing import Optional, Dict, List

import yaml
import pandas as pd
from modules.atomistic.gromacs.parser.itp_parser import ITPParser
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from rdkit_new.base_polymer_generator import BasePolymerGenerator
from typing import Optional, Dict, List


class OpenMSCGMappingGenerator:
    """
    Generates an OpenMSCG-compatible YAML mapping file.
    """

    def __init__(
        self,
        polymer: BasePolymerGenerator,
        gro_file: str,
        solvent_itp: Optional[str] = None,
        molecule_name: Optional[str] = "UNL",
    ):
        """
        Initializes the OpenMSCG YAML generator.

        :param polymer: Polymer generator object containing CG bead mappings.
        :param gro_file: Path to the GROMACS .gro file (needed for anchor alignment).
        :param solvent_itp: Path to the solvent .itp file (optional).
        :param molecule_name: Name of the molecule (default: "UNL").
        """
        self.polymer = polymer
        self.gro_file = gro_file
        self.molecule_name = molecule_name
        self.bead_mappings = polymer.cg_map
        self.solvent_mapping = (
            self._parse_solvent_itp(solvent_itp) if solvent_itp else {}
        )
        self.residue_map = self._parse_gro_residues(gro_file)

    def _parse_gro_residues(self, gro_file: str) -> Dict[int, int]:
        """
        Parses the .gro file to extract residue indices.

        :param gro_file: Path to the .gro file.
        :return: Dictionary mapping residue number to atom indices.
        """
        sections = GromacsParser().parse(gro_file)
        gro_section = next(iter(sections.values()))
        gro_handler = GroHandler()
        gro_handler.process(gro_section)
        gro_df = gro_handler.content

        # Extract unique residue numbers and their first atom index
        residue_map = {}
        for _, row in gro_df.iterrows():
            res_number = int(row["Residue Number"])  # Residue number
            atom_index = int(row["Atom Index"]) - 1  # Convert 1-based to 0-based
            if res_number not in residue_map:
                residue_map[res_number] = atom_index  # Store first occurrence

        return residue_map

    def _parse_solvent_itp(self, itp_file: str) -> Dict:
        """
        Parses the solvent .itp file to extract residue and atom names.

        :param itp_file: Path to the solvent .itp file.
        :return: Dictionary with solvent bead mappings.
        """
        parser = ITPParser(itp_file)
        itp_df = parser.retrieve_content("atoms")

        # Extract the solvent residue name
        solvent_res = itp_df["res"].iloc[0]

        # Map all solvent atoms to a single bead
        return {
            "name": solvent_res,
            "index": list(range(len(itp_df))),  # Indexing all atoms in solvent
            "x-weight": [1.0] * len(itp_df),  # Equal weight distribution
            "f-weight": [1.0] * len(itp_df),
        }

    def generate_mapping_yaml(self, output_file: str):
        """
        Generates the YAML mapping file.

        :param output_file: Output file path.
        """
        mapping_data = {"site-types": {}, "system": []}

        # === ADD POLYMER BEAD MAPPINGS ===
        monomer_size = max(len(bead["atom_indices"]) for bead in self.bead_mappings)
        polymer_residues = sorted(self.residue_map.keys())  # Sorted residue numbers
        polymer_repeat = len(polymer_residues)

        for bead in self.bead_mappings:
            bead_name = bead["unique_name"]
            atom_indices = list(bead["atom_indices"])
            x_weights = list(bead["x-weight"])
            f_weights = list(bead["f-weight"])

            mapping_data["site-types"][bead_name] = {
                "index": atom_indices,
                "x-weight": x_weights,
                "f-weight": f_weights,
            }

        # Determine polymer anchor (first residue in .gro)
        polymer_anchor = self.residue_map[polymer_residues[0]]

        # Define polymer system mapping
        polymer_entry = {
            "anchor": polymer_anchor,
            "repeat": polymer_repeat,  # Number of polymer units
            "offset": monomer_size,  # Atoms per monomer
            "sites": [
                [bead["unique_name"], i * monomer_size]
                for i, bead in enumerate(self.bead_mappings)
            ],
        }
        mapping_data["system"].append(polymer_entry)

        # === ADD SOLVENT MAPPING (IF EXISTS) ===
        if self.solvent_mapping:
            solvent_anchor = self.residue_map[
                max(polymer_residues) + 1
            ]  # First solvent residue
            solvent_entry = {
                "anchor": solvent_anchor,  # Solvent starts after polymer atoms
                "repeat": 1,  # One bead per solvent molecule
                "offset": len(self.solvent_mapping["index"]),
                "sites": [[self.solvent_mapping["name"], 0]],
            }
            mapping_data["system"].append(solvent_entry)

        # === WRITE YAML WITH FIXED FORMATTING ===
        with open(output_file, "w") as f:
            yaml.dump(
                mapping_data,
                f,
                default_flow_style=False,
                sort_keys=False,
                allow_unicode=True,
            )

        print(f"OpenMSCG mapping file saved to {output_file}")


# === USAGE EXAMPLE ===
from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator

polymer = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
polymer.generate_polymer(3, "rdkit_test2", overwrite=False, save=False)

# Example Usage:
generator = OpenMSCGMappingGenerator(
    polymer, gro_file="temp/production.gro", solvent_itp="1_5_test/hexane/solvent.itp"
)
generator.generate_mapping_yaml("cgmap.yaml")
