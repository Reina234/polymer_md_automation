import os
from config.gromacs_paths import get_template_path, get_script_path

class TemplateProcessor:
    """
    Processes GROMACS templates by replacing placeholders with specific values.
    """

    @staticmethod
    def process_template(step: str, temperature: int, replacements: dict):
        """
        Process a base template for a specific GROMACS step and temperature by replacing placeholders,
        and save the processed script in the appropriate directory.

        Args:
            step (str): Simulation step (e.g., "em", "nvt", "npt").
            temperature (int): Temperature in Kelvin.
            replacements (dict): Dictionary of placeholder keys and their replacements.

        Returns:
            str: Path to the processed script file.
        """
        # Get the base template path
        base_template_path = get_template_path(step)

        # Get the path where the temperature-specific script will be saved
        script_path = get_script_path(step, temperature)

        if not os.path.exists(base_template_path):
            raise FileNotFoundError(f"Base template file not found: {base_template_path}")

        with open(base_template_path, "r") as file:
            content = file.read()

        # Replace placeholders
        for key, value in replacements.items():
            content = content.replace(f"@@{key}@@", str(value))

        # Ensure the directory exists
        os.makedirs(os.path.dirname(script_path), exist_ok=True)

        # Write the processed content to the script file
        with open(script_path, "w") as file:
            file.write(content)

        return script_path
