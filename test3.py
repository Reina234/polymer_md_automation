from modules.workflows.separated.parametiser.solvent import SolventParametiser


generator = SolventParametiser(
    solvent_name="acetic_acid",
    solvent_smiles="CC(=O)O",
    solvent_compressibility="0.00045",
    solvent_density=1.049,
    output_dir="test",
).run()

from modules.workflows.separated.parametiser.polymer import PolymerParametiser

generator = PolymerParametiser(monomer_smiles=["C=C"], num_units=3, output_dir="test")
generator.run()

from zzz_directory_parser.polymer_directory_parser import PolymerDirectoryParser
from zzz_directory_parser.solvent_directory_parser import SolventDirectoryParser

polymer_parser = PolymerDirectoryParser("test").parse_directory()
solvent_parser = SolventDirectoryParser("test").parse_directory()
print(polymer_parser)
print(solvent_parser)

parameterised_solvent, solvent, solvent_smiles = solvent_parser[0]
parameterised_polymer, monomer_smiles, n_units = polymer_parser[0]

from modules.workflows.separated.gromacs.joined import JoinedAtomisticPolymerWorkflow

workflow = JoinedAtomisticPolymerWorkflow(
    parameterised_polymer=parameterised_polymer,
    monomer_smiles=monomer_smiles,
    num_units=n_units,
    solvent=solvent,
    solvent_smiles=solvent_smiles,
    parameterised_solvent=parameterised_solvent,
    temperature=298,
    output_dir="test3",
    csv_file_path="output_test",
).run()
