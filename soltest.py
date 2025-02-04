from modules.rdkit.solvent_generator import SolventGenerator

from modules.workflows.atomistic.joined_workflow import JoinedAtomisticPolymerWorkflow


workflow = JoinedAtomisticPolymerWorkflow(
    monomer_smiles=["C=Cc1ccccc1", "C=C"],
    num_units=10,
    solvent_name="hexane",
    solvent_smiles="CCCCCC",
    solvent_density=661,
    temperatures=[298],
    output_dir="ZZZ",
    solvent_compressibility=1.24e-4,
    csv_file_path="test",
)

workflow.run()
