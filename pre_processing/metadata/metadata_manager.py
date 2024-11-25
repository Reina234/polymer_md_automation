import os
class MetadataManager:
    """
    Tracks and manages metadata for the workflow.
    """

    def __init__(self):
        self.metadata = {"components": []}
        self.notes = ""

    def register_component(self, name: str, version: str, description: str, options: dict):
        """
        Register metadata for a component.
        """
        self.metadata["components"].append({
            "name": name,
            "version": version,
            "description": description,
            "options": options
        })

    def add_notes(self, notes: str):
        """
        Add custom notes to the metadata.
        """
        self.notes = notes

    def save_metadata(self, output_dir: str):
        """
        Save metadata to a file.
        """
        metadata_file = os.path.join(output_dir, "metadata.json")
        metadata_to_save = {
            "workflow_metadata": self.metadata,
            "user_notes": self.notes,
        }
        with open(metadata_file, "w") as f:
            import json
            json.dump(metadata_to_save, f, indent=4)
