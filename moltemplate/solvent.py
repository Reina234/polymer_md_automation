from moltemplate.base_molecule import BaseMoltemplateMolecule
from typing import Dict


class MoltemplateSolvent(BaseMoltemplateMolecule):
    def __init__(self, open_mscg_mass_map=Dict[str, float], sol_resname="SOL"):
        super().__init__(open_mscg_mass_map=open_mscg_mass_map)
        self.resname = sol_resname
        self.mass = self.premade_mass_map[self.resname]
        self.mass_mapping[self.resname] = self.mass

    def generate_lt(self) -> str:
        return f"Molecule Solvent {{\n  $atom:{self.resname} 0.0 0.0 0.0\n}}\n"
