# File: processing/common/metadata_tracker.py

import os
import json
from A_modules.shared.utils.file_utils import check_directory_exists
import logging

logger = logging.getLogger(__name__)


class MetadataTracker:
    """
    Tracks and manages metadata for the workflow.
    """

    FILE_EXTENSION = "json"
    FILE_NAME = "metadata"

    def __init__(self):
        self.metadata = {"steps": []}

    def add_step(self, step_name: str, details: dict):
        """
        Add a step to the metadata.

        :param step_name: Name of the step
        :type step_name: str
        :param details: Details of the step
        :type details: dict
        """
        self.metadata["steps"].append({"name": step_name, "details": details})

    def save(
        self,
        output_dir: str,
        file_name: str = FILE_NAME,
        file_extension: str = FILE_EXTENSION,
    ):
        output_dir = check_directory_exists(output_dir)
        metadata_file_path = os.path.join(output_dir, f"{file_name}.{file_extension}")

        with open(metadata_file_path, "w") as f:
            json.dump(self.metadata, f, indent=4)

        logging.info(f"Metadata saved to: {metadata_file_path}")
