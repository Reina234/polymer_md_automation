from input_data.monomer_smiles import monomer_smiles_list
from separated_simulation_manager import PolymerParametiser
import os

job_id = os.getenv("SLURM_JOB_ID", "local_run")

if __name__ == "__main__":
    manager = PolymerParametiser(
        monomer_smiles=monomer_smiles_list,
        num_units=[5, 10, 20],
        output_dir=f"outputs_{job_id}",
    )

    manager.run()
