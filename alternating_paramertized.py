from rdkit_new.homopolymer_generator import HomopolymerGenerator
from config.paths import TEMP_DIR
from config.acpype_config import AcpypeOutputConfig
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from modules.shared.file_conversion.converter_factory import ConverterFactory
from modules.shared.utils.file_utils import copy_file, delete_directory
from modules.atomistic.acpype_parameterizer.acpype_parametizer import (
    ACPYPEParameterizer,
)
from data_models.output_types import GromacsPaths
from config.paths import PARAMETERISED_POLYMER_DIR
import os

# NOTE: I think I can use just a copolymer generator, and create a homopolymer by passing in two of the same monomer units


def parameterize(
    polymer_pdb: str,
    output_dir: str = PARAMETERISED_POLYMER_DIR,
    temp_dir: str = TEMP_DIR,
    keep_pdb: bool = False,
    keep_mol2: bool = False,
    polymer_name: str = "POLY",
    parameterizer: ACPYPEParameterizer = ACPYPEParameterizer,
    converter_factory: ConverterFactory = ConverterFactory(),
    cleanup: bool = True,
) -> GromacsPaths:

    if keep_pdb:
        polymer_pdb = copy_file(polymer_pdb, output_dir)

    if keep_mol2:
        mol2_output_dir = output_dir
    else:
        mol2_output_dir = temp_dir
    mol2_converter = converter_factory.get_converter("pdb", "mol2")
    mol2_file = mol2_converter.run(polymer_pdb, mol2_output_dir)
    parameterizer = parameterizer(acpype_molecule_name=polymer_name)
    file_config = AcpypeOutputConfig(itp=True, gro=True, top=True, posre=False)
    parametized_polymer = parameterizer.run(mol2_file, output_dir, file_config)
    if cleanup:
        delete_directory(temp_dir, confirm=False)
    return parametized_polymer


from rdkit_new.alternating_copolymer import AlternatingPolymerGenerator
from rdkit import Chem

from rdkit_new.votca_mapping_generator import VOTCAMappingGenerator
from modules.atomistic.gromacs.parser.gromacs_parser import GromacsParser
from modules.atomistic.gromacs.parser.handlers.data_handler import DataHandler
from rdkit import Chem
import pandas as pd
from typing import Dict, List, Tuple

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator.generate_polymer(10, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

import numpy as np
import MDAnalysis as mda
from collections import defaultdict


import pandas as pd
import re


class PolymerITPExtender:
    num_cols = {
        "bonds": ["ai", "aj"],
        "pairs": ["ai", "aj"],
        "angles": ["ai", "aj", "ak"],
        "dihedrals": ["ai", "aj", "ak", "al"],
        "impropers": ["ai", "aj", "ak", "al"],
        "atoms": ["nr", "cgnr"],
    }
    SECTION_HEADERS = {
        "atoms": ["nr", "type", "resnr", "residue", "atom", "cgnr", "charge", "mass"],
        "bonds": ["ai", "aj", "funct", "r", "k"],
        "pairs": ["ai", "aj", "funct"],
        "angles": ["ai", "aj", "ak", "funct", "theta", "k_theta"],
        "dihedrals": ["ai", "aj", "ak", "al", "funct", "phi", "k_phi", "pn"],
        "impropers": ["ai", "aj", "ak", "al", "funct", "phi", "k_phi", "pn"],
    }
    comment_sections = {
        k: v for k, v in num_cols.items() if k != "atoms"
    }  # Exclude "atoms"

    def __init__(self, itp_path, cg_map, n_repeat, atom_start_index=None):
        """
        Args:
            itp_path (str): Path to .itp file.
            cg_map (list): Mapping of residue indices.
            n_repeat (int): Number of times to repeat the middle section.
            atom_start_index (int or None): Determines atom name counting style.
                - None → "C", "H", etc.
                - 0 → "C0", "H0", etc.
                - 1 → "C1", "H1", etc.
        """
        self.itp_path = itp_path
        self.cg_map = cg_map
        self.n_repeat = n_repeat
        self.atom_start_index = atom_start_index  # Determines how atom indices start
        self.sections = {}

        self._middle_length = self._calculate_middle_length()
        self._parse_itp()
        self._process_sections()

    def _parse_itp(self) -> None:
        """Parses .itp file into DataFrames while preserving interaction types."""
        itp_file = self.itp_path
        sections = {}
        current_section = None
        section_data = []

        with open(itp_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                # **Remove in-line comments** before processing
                line = re.split(r"\s*;\s*", line)[0]

                # Detect section headers
                match = re.match(r"^\[\s*(\w+)\s*\](.*)$", line)
                if match:
                    if current_section and section_data:
                        sections[current_section] = self._create_dataframe(
                            current_section, section_data
                        )
                        section_data = []

                    section_name, extra_info = match.groups()
                    section_name = section_name.lower()

                    if "improper" in extra_info.lower():
                        current_section = "impropers"
                    else:
                        current_section = section_name

                    continue

                if current_section:
                    section_data.append(line.split())

        # Convert last section into DataFrame
        if current_section and section_data:
            sections[current_section] = self._create_dataframe(
                current_section, section_data
            )

        self.sections = sections

    def _create_dataframe(self, section_name, data):
        """Creates a DataFrame with predefined headers, handling extra/missing columns."""
        headers = self.SECTION_HEADERS.get(section_name, None)
        if headers:

            df = pd.DataFrame(data, columns=headers)

            # Fill missing headers with NaN
            for col in self.SECTION_HEADERS.get(section_name, []):
                if col not in df.columns:
                    df[col] = None
        else:
            df = pd.DataFrame(data)  # Generic DataFrame for unknown sections

        return df

    def _calculate_middle_length(self):
        """Determines the number of atoms in the repeating middle section."""
        middle_atoms = sum(len(entry["atom_indices"]) for entry in self.cg_map[1:-1])
        return middle_atoms

    def _shift_atom_indices(self):
        """Shifts atom names based on the defined counting style."""
        if "atoms" not in self.sections:
            return

        df = self.sections["atoms"]
        atom_counts = {}  # Track occurrences of each atom type

        for i, row in df.iterrows():
            atom_name = row["atom"]
            match = re.match(
                r"([A-Za-z]+)(\d*)", atom_name
            )  # Extract base name & number
            if match:
                base_name, number = match.groups()
                if base_name not in atom_counts:
                    atom_counts[base_name] = self.atom_start_index or 0  # Start count

                df.at[i, "atom"] = f"{base_name}{atom_counts[base_name]}"
                atom_counts[base_name] += 1  # Increment count

        self.sections["atoms"] = df

    def _process_sections(self):
        """Processes all relevant sections, shifts indices, and duplicates the middle section."""
        num_cols = {
            "bonds": ["ai", "aj"],
            "pairs": ["ai", "aj"],
            "angles": ["ai", "aj", "ak"],
            "dihedrals": ["ai", "aj", "ak", "al"],
            "impropers": ["ai", "aj", "ak", "al"],
            "atoms": ["nr", "cgnr"],
        }

        for section, cols in num_cols.items():
            if section in self.sections:
                df = self.sections[section]

                # Ensure numeric conversion
                for col in cols:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

                # Identify start, middle, and end
                start_indices = set(self.cg_map[0]["atom_indices"])
                middle_indices = set(
                    idx for entry in self.cg_map[1:-1] for idx in entry["atom_indices"]
                )

                start_df = df[df[cols[0]].isin(start_indices)].copy()
                middle_df = df[df[cols[0]].isin(middle_indices)].copy()

                assigned_indices = start_indices.union(middle_indices)
                all_indices = set(df[cols[0]].unique())

                # End should contain everything that isn't in start or middle
                end_indices = all_indices - assigned_indices
                end_df = df[df[cols[0]].isin(end_indices)].copy()

                # Duplicate middle n times and shift indices accordingly
                all_middle = []
                for i in range(self.n_repeat + 1):
                    temp_df = middle_df.copy()
                    for col in cols:
                        temp_df[col] += (
                            i * self._middle_length
                        )  # Shift middle section indices
                    all_middle.append(temp_df)

                # 🔄 Apply index shift to `end_df`
                for col in cols:
                    end_df[col] += (
                        self.n_repeat * self._middle_length
                    )  # Shift end section indices

                # Combine everything back together
                df = pd.concat([start_df] + all_middle + [end_df], ignore_index=True)
                self.sections[section] = df

        self._shift_atom_indices()
        self._add_comments()

    def _add_comments(self):
        """Adds comments to relevant sections based on atom names."""
        if "atoms" not in self.sections:
            return

        atom_map = {
            row["nr"]: row["atom"] for _, row in self.sections["atoms"].iterrows()
        }  # nr -> atom_name

        for section, cols in self.comment_sections.items():
            if section in self.sections:
                df = self.sections[section]
                df["comment"] = df.apply(
                    lambda row: ";"
                    + " - ".join(atom_map.get(row[col], "?") for col in cols),
                    axis=1,
                )

    def write_itp(self, output_path):
        """Writes the modified ITP file to the specified output path."""
        with open(output_path, "w") as f:
            for section, df in self.sections.items():
                f.write(f"\n[ {section} ]\n")
                for _, row in df.iterrows():
                    line = " ".join(map(str, row.dropna().tolist()))
                    f.write(f"{line}\n")


generator2 = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator2.generate_polymer(3, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

itp = PolymerITPExtender(
    "1_22_test/c=cc1ccccc1_3/POLY_GMX.itp",
    generator2.cg_map,
    7,
)

itp.write_itp("test_3.itp")
