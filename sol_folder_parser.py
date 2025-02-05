import os
from dataclasses import dataclass, field
from typing import Optional, List, Tuple, Dict


class SolventDirParser:
    acpype_output_basename = "SVT_GMX"

    def __init__(self, input_dir: str):
        self.input_dir = input_dir

    @staticmethod
    def parse_solvent_properties(file_path: str) -> Dict[str, str]:
        properties = {}

        if not os.path.exists(file_path):
            return properties  # Return empty dictionary if file doesn't exist

        with open(file_path, "r") as file:
            for line in file:
                line = line.strip()
                if "=" in line:
                    key, value = line.split("=", 1)  # Only split on first '='
                    properties[key.strip()] = value.strip()

        return properties


def parse_directory(base_dir: str) -> List[Tuple[GromacsPaths, Solvent]]:
    """
    Iterates through the base directory to extract solvent properties and GROMACS file paths.

    :param base_dir: The base directory containing solvent subdirectories.
    :return: A list of tuples (GromacsPaths instance, Solvent instance).
    """
    results = []

    for subdir in os.listdir(base_dir):
        subdir_path = os.path.join(base_dir, subdir)
        if not os.path.isdir(subdir_path):
            continue  # Skip files, only process directories

        # Define expected file paths
        itp_file = os.path.join(subdir_path, "SVT_GMX.itp")
        gro_file = os.path.join(subdir_path, "SVT_GMX.gro")
        top_file = os.path.join(subdir_path, "SVT_GMX.top")
        sol_properties_file = os.path.join(subdir_path, "sol_properties.txt")

        # Extract solvent properties
        properties = parse_solvent_properties(sol_properties_file)

        if not properties:
            continue  # Skip if solvent properties file is missing or empty

        # Create GromacsPaths instance
        gromacs_paths = GromacsPaths(
            itp_path=itp_file if os.path.exists(itp_file) else None,
            gro_path=gro_file if os.path.exists(gro_file) else None,
            top_path=top_file if os.path.exists(top_file) else None,
        )

        # Ensure required fields exist
        required_fields = {"name", "compressibility", "density", "pdb_path", "SMILES"}
        if not required_fields.issubset(properties.keys()):
            print(
                f"Warning: Missing required fields in {sol_properties_file}. Skipping..."
            )
            continue

        # Convert numeric fields safely
        try:
            compressibility = float(properties["compressibility"])
            density = float(properties["density"])
        except ValueError:
            print(
                f"Error: Invalid numeric values in {sol_properties_file}. Skipping..."
            )
            continue

        # Create Solvent instance
        solvent = Solvent(
            name=properties["name"],
            compressibility=compressibility,
            density=density,
            pdb_path=properties["pdb_path"],
            SMILES=properties["SMILES"],
            molecular_weight=(
                float(properties["molecular_weight"])
                if "molecular_weight" in properties
                else None
            ),
        )

        # Append result tuple
        results.append((gromacs_paths, solvent))

    return results
