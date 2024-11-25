import os

class MetadataTracker:
    """
    Tracks metadata across the workflow, ensuring propagation between processes.
    """

    def __init__(self):
        self.metadata = {
            "processes": [],
            "history": []
        }

    def add_process_metadata(self, name: str, version: str, description: str, options: dict):
        """
        Add metadata for a process/tool.
        """
        process_metadata = {
            "name": name,
            "version": version,
            "description": description,
            "options": options
        }
        self.metadata["processes"].append(process_metadata)
        self.metadata["history"].append(name)

    def get_metadata(self):
        """
        Retrieve the complete metadata.
        """
        return self.metadata

    def save_metadata(self, output_dir: str):
        """
        Save metadata to a file.
        """
        import json
        metadata_file = os.path.join(output_dir, "metadata.json")
        with open(metadata_file, "w") as f:
            json.dump(self.metadata, f, indent=4)
