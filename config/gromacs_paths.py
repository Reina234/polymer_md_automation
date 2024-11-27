import os

# Base directories
BASE_DIR = os.path.abspath(os.path.dirname(__file__))  # Adjust based on your project structure
TEMP_DIR = os.path.join(BASE_DIR, "temp")
SCRIPTS_DIR = os.path.join(BASE_DIR, "scripts")

GROMACS_SCRIPTS_DIR = os.path.join(SCRIPTS_DIR, "gromacs")
GROMACS_TEMPLATES_DIR = os.path.join(GROMACS_SCRIPTS_DIR, "templates")

# Base templates with no temperature placeholders
GROMACS_TEMPLATES = {
    "ions": os.path.join(GROMACS_TEMPLATES_DIR, "ions.mdp"),
    "em": os.path.join(GROMACS_TEMPLATES_DIR, "em.mdp"),
    "nvt": os.path.join(GROMACS_TEMPLATES_DIR, "nvt.mdp"),
    "npt": os.path.join(GROMACS_TEMPLATES_DIR, "npt.mdp"),
}

################################## NOTE: NEEDS CHANGING ##################################
GROMACS_ION_SCRIPT = GROMACS_TEMPLATES["ions"]
##########################################################################################


# Directories where temperature-specific scripts are located
GROMACS_EM_DIR = os.path.join(GROMACS_SCRIPTS_DIR, "em")
GROMACS_NVT_DIR = os.path.join(GROMACS_SCRIPTS_DIR, "nvt")
GROMACS_NPT_DIR = os.path.join(GROMACS_SCRIPTS_DIR, "npt")

# Paths to temperature-specified scripts with placeholders
SCRIPT_PATHS = {
    "em": os.path.join(GROMACS_EM_DIR, "em_@@TEMPERATURE@@K.mdp"),
    "nvt": os.path.join(GROMACS_NVT_DIR, "nvt_@@TEMPERATURE@@K.mdp"),
    "npt": os.path.join(GROMACS_NPT_DIR, "npt_@@TEMPERATURE@@K.mdp"),
}

# Other fixed paths
SOLVENTS_DIR = os.path.join(BASE_DIR, "solvents")
LOGS_DIR = os.path.join(BASE_DIR, "logs")


def get_template_path(step: str) -> str:
    """
    Get the path to the base template for a specific step.

    Args:
        step (str): The simulation step (e.g., "em", "nvt", "npt", "ions").

    Returns:
        str: The path to the base template file.
    """
    if step not in GROMACS_TEMPLATES:
        raise ValueError(f"Unknown simulation step: {step}. Valid options: {list(GROMACS_TEMPLATES.keys())}")
    return GROMACS_TEMPLATES[step]


def get_script_path(step: str, temperature: int) -> str:
    """
    Get the path to the temperature-specific script for a specific step.

    Args:
        step (str): The simulation step (e.g., "em", "nvt", "npt").
        temperature (int): The temperature in Kelvin.

    Returns:
        str: The resolved script path.
    """
    if step not in SCRIPT_PATHS:
        raise ValueError(f"Unknown simulation step: {step}. Valid options: {list(SCRIPT_PATHS.keys())}")
    return SCRIPT_PATHS[step].replace("@@TEMPERATURE@@", str(temperature))


def get_output_dir(step: str, temperature: int) -> str:
    """
    Get the output directory for a specific step and temperature.

    Args:
        step (str): The simulation step (e.g., "em", "nvt", "npt", "ions").
        temperature (int): The temperature in Kelvin.

    Returns:
        str: The resolved output directory path.
    """
    return os.path.join(TEMP_DIR, step, f"{temperature}K")

