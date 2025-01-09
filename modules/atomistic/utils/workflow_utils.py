import re
import os
from typing import List, Dict
from config.paths import EQUILIBRIATED_SOLVENT_BOX_DIR, EQUILIBRIATED_OUTPUTS_SUBDIR

import os


def iterate_files_with_suffix(directory, suffix):
    """
    Iterate through files in a directory with a specific suffix.

    :param directory: The path to the directory.
    :param suffix: The file suffix to filter (e.g., ".gro").
    :return: A generator of file paths.
    """
    for file_name in os.listdir(directory):
        if file_name.endswith(suffix):
            yield os.path.join(directory, file_name)


def retrieve_temp_and_compressibility_from_dir(
    solvent_name: str, solvent_box_dir: str = EQUILIBRIATED_SOLVENT_BOX_DIR
) -> Dict[str, str]:
    subdir = EQUILIBRIATED_OUTPUTS_SUBDIR
    gro_dir = os.path.join(solvent_box_dir, solvent_name, subdir)
    list_of_params = []
    for file_path in iterate_files_with_suffix(gro_dir, ".gro"):
        list_of_params.append(
            extract_specific_params(file_path, ["temp", "compressibility"])
        )
    return list_of_params


def extract_specific_params(
    filename: str, params_to_extract: List[str] = ["temp", "compressibility"]
) -> Dict[str, str]:
    """
    Extracts specific parameter values from a filename.

    :param filename: The filename in the format 'param1_value1_param2_value2.extension'
    :param params_to_extract: A list of parameter names to extract (e.g., ['temp', 'compressibility']).
    :return: A dictionary of extracted parameter values.
    """
    # Create a regex pattern for extracting key-value pairs
    pattern = r"(?P<key>\w+)_(?P<value>[0-9.eE-]+)"

    matches = re.finditer(pattern, filename)
    extracted_params = {}

    for match in matches:
        key = match.group("key")
        value = match.group("value")
        if key in params_to_extract:
            extracted_params[key] = value

    return extracted_params
