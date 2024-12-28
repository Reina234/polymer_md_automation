from dataclasses import dataclass
from typing import Optional, List


@dataclass
class GromacsPaths:
    """
    A class to store the paths of the acpype files.
    """

    itp_path: Optional[str] = None
    gro_path: Optional[str] = None
    top_path: Optional[str] = None
    posre_path: Optional[str] = None

    def to_list(self) -> List[Optional[str]]:
        """
        Converts the GromacsPaths instance into a list of file paths.

        :return: A list of file paths.
        :rtype: List[Optional[str]]
        """
        return [self.itp_path, self.gro_path, self.top_path, self.posre_path]
