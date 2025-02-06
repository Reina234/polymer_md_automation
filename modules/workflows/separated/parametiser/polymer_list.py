from modules.workflows.separated.parametiser.polymer import PolymerParametiser
import itertools
from typing import List


class PolymerListParametiser:
    def __init__(
        self, full_smiles_list: list, output_dir: str, num_units_list: List[int]
    ):
        self.full_smiles_list = full_smiles_list
        self.output_dir = output_dir
        self.num_units_list = num_units_list
        self.monomer_list_combinations = self._generate_combinations()

    def _generate_combinations(self):
        all_combinations = []
        for monomer in self.full_smiles_list:
            all_combinations.append([monomer])

        for monomer_pair in itertools.combinations(self.full_smiles_list, 2):
            all_combinations.append(list(monomer_pair))

        for monomer_triplet in itertools.combinations(self.full_smiles_list, 3):
            all_combinations.append(list(monomer_triplet))

        return all_combinations

    def run(self):
        for num_units in self.num_units_list:
            for monomer_list in self.monomer_list_combinations:
                generator = PolymerParametiser(
                    monomer_smiles=monomer_list,
                    num_units=num_units,
                    output_dir=self.output_dir,
                )
                generator.run()
