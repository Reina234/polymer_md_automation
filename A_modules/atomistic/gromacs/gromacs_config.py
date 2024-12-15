from enum import Enum

MDP_NAMING_SCHEME = "{run_type}_{temp}K.mdp"
BASE_TEMPLATE_PATH = "A_modules/atomistic/gromacs/mdp_templates"


class TemplatedMdps(Enum):
    NVT = ("nvt", f"{BASE_TEMPLATE_PATH}/nvt.mdp")
    NPT = ("npt", f"{BASE_TEMPLATE_PATH}/npt.mdp")
    MINIM = ("minim", f"{BASE_TEMPLATE_PATH}/minim.mdp")

    def __init__(self, run_type: str, path: str):
        self.run_type = run_type  # Avoid using `name` to prevent conflict
        self.path = path

    @property
    def template_path(self) -> str:
        return self.path

    def generate_mdp_name(
        self, temp: int, mdp_output_naming_scheme: str = MDP_NAMING_SCHEME
    ) -> str:
        """
        Generate an MDP file name based on the naming scheme and temperature.

        :param temp: Temperature in Kelvin to include in the file name.
        :type temp: int
        :return: Generated MDP file name.
        :rtype: str
        """
        return mdp_output_naming_scheme.format(run_type=self.run_type, temp=temp)


# NOTE: will need to change this to the correct path after renaming
