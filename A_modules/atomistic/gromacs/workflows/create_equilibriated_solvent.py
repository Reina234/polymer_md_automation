from A_modules.atomistic.gromacs.commands.base_gromacs_command import BaseGromacsCommand
from A_config.constants import MassUnits2, LengthUnits2
from A_modules.atomistic.gromacs.commands.editconf import Editconf
from A_modules.atomistic.gromacs.commands.solvate import Solvate
from A_modules.shared.utils.calculation_utils import calculate_num_particles
from A_modules.shared.metadata_tracker import MetadataTracker
from A_modules.atomistic.gromacs.utils.utils import calculate_minimum_box_size_from_gro
from typing import Optional, List, Dict

output_names = {
    "editconf_solvent_box": "editconf_solvent_box.gro",
    "editconf_box_dim": "editconf_box_dim.gro",
    "solvated_box": "solvent_box.gro",
}


def test(
    gro_file: str,
    output_dir: str,
    input_topol_path: str,
    density: float,
    box_size: List[float],
    molecular_weight: float,
    names: Dict[str, str] = output_names,
    metadata_tracker: Optional[MetadataTracker] = None,
    box_units: LengthUnits2 = LengthUnits2.NANOMETER,
    mw_units: MassUnits2 = MassUnits2.GRAM,
    editconf: Editconf = Editconf(),
    solvate: Solvate = Solvate(),
    padding: float = 0.1,
    initial_box_factor: float = 0.9,
):
    # NOTE: will need a conversion method if the units are not nw

    # ADD IN  the solvated box
    editconf_solvent_box = editconf.run(
        gro_file,
        output_dir,
        box_size_nm=initial_box_size,
        output_name=output_names["editconf_solvent_box"],
    )

    solvated_box = solvate.run(
        editconf_solvent_box,
        editconf_box_dim,
        input_topol_path=input_topol_path,
        output_dir=output_dir,
        output_name=output_names["solvated_box"],
    )


# define loop
