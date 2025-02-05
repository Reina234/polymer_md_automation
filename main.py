from input_data.monomer_smiles import monomer_smiles_list
from simulation_manager import PolymerSimulationManager
import os

job_id = os.getenv("SLURM_JOB_ID", "local_run")

if __name__ == "__main__":
    manager = PolymerSimulationManager(
        solvent_csv="input_data/solvent_data.csv",
        monomer_smiles=monomer_smiles_list,
        num_units=[5, 10, 20],
        temperatures=[280, 298, 348],
        output_dir=f"outputs_{job_id}",
        csv_file_path="outputs",
    )

    manager.run()
