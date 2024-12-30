from typing import Dict, List
from modules.shared.packmol.base_packmol_operation import BasePackmolOperation
from data_models.solvent import Solvent
from modules.shared.utils.calculation_utils import (
    calculate_num_particles,
    box_dimensions_check_wrapper,
)
from config.constants import LengthUnits, MassUnits
from modules.shared.utils.file_utils import file_type_check_wrapper


class PackmolSolventBox(BasePackmolOperation):
    template_name = "solvent_box_template.inp"

    def generate_template_params(
        self,
        input_pdb: str,
        output_file: str,
        solvent: Solvent,
        box_size_nm: List[float],
        **kwargs,
    ) -> Dict[str, str]:

        box_size_angstrom = [
            dim * LengthUnits.NANOMETER.value / LengthUnits.ANGSTROM.value
            for dim in box_size_nm
        ]

        print("box_size_angstrom", box_size_angstrom)
        num_molecules = calculate_num_particles(
            box_dimensions=box_size_angstrom,
            molecular_weight=solvent.molecular_weight,
            density_SI=solvent.density,
            box_units=LengthUnits.ANGSTROM,
            mass_units=MassUnits.GRAM,
        )
        print(num_molecules)
        return {
            "output_file": output_file,
            "solvent_file": input_pdb,
            "num_molecules": num_molecules,
            "box_size_x": box_size_angstrom[0],
            "box_size_y": box_size_angstrom[1],
            "box_size_z": box_size_angstrom[2],
        }
