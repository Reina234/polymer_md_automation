from dataclasses import dataclass
from typing import Optional, List


@dataclass
class AcpypePaths:
    """
    A class to store the paths of the acpype files.
    """

    itp_path: Optional[str] = None
    gro_path: Optional[str] = None
    top_path: Optional[str] = None
    posre_path: Optional[str] = None

    def to_list(self) -> List[Optional[str]]:
        """
        Converts the AcpypePaths instance into a list of file paths.

        :return: A list of file paths.
        :rtype: List[Optional[str]]
        """
        return [self.itp_path, self.gro_path, self.top_path, self.posre_path]


@dataclass
class AcpypeOutputConfig:
    """
    A class to store the configuration for the acpype output.
    """

    itp: bool
    gro: bool
    top: bool
    posre: bool

    ITP_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.itp"
    GRO_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.gro"
    TOP_FILE_NAME_FORMAT: str = "{molecule_name}_GMX.top"
    POSRE_FILE_NAME_FORMAT: str = "posre_{molecule_name}.itp"
