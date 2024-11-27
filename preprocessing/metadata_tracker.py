# File: processing/common/metadata_tracker.py

import os
import json


class MetadataTracker:
    """
    Tracks and manages metadata for the workflow.
    """

    def __init__(self):
        self.metadata = {"steps": []}

    def add_step(self, step_name: str, details: dict):
        """
        Add a step to the metadata.

        Args:
            step_name (str): Name of the step (e.g., "ACPYPE Parameterization").
            details (dict): Details about the step (e.g., input/output files, parameters).
        """
        self.metadata["steps"].append({"name": step_name, "details": details})

    def save(self, output_dir: str, filename: str = "metadata.json"):
        """
        Save metadata to a file in the specified output directory.

        Args:
            output_dir (str): Directory where the metadata file will be saved.
            filename (str): Name of the metadata file (default: "metadata.json").
        """
        # Ensure the output directory exists
        os.makedirs(output_dir, exist_ok=True)

        # Construct the full path for the metadata file
        metadata_file_path = os.path.join(output_dir, filename)

        # Save metadata as JSON
        with open(metadata_file_path, "w") as f:
            json.dump(self.metadata, f, indent=4)

        print(f"Metadata saved to: {metadata_file_path}")
