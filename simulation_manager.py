import itertools
import logging
import os
import pandas as pd
from typing import List, Tuple
from modules.workflows.atomistic.joined_workflow import JoinedAtomisticPolymerWorkflow


# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PolymerSimulationManager:
    def __init__(
        self,
        solvent_csv: str,
        monomer_smiles: List[str],
        output_dir: str,
        progress_file: str = "progress",
        csv_file_path: str = "outputs",
        num_units: List[int] = [5, 10, 20],
        temperatures: List[int] = [280, 298, 346],
    ):
        """
        Initializes the simulation manager.
        :param solvent_csv: Path to the CSV file containing solvent properties.
        :param monomer_smiles: List of monomer SMILES strings.
        :param output_dir: Directory where results should be stored.
        :param progress_file: File for tracking progress.
        """
        self.solvent_df = pd.read_csv(solvent_csv)
        self.monomer_smiles = monomer_smiles
        self.num_units = num_units
        self.temperatures = temperatures
        self.output_dir = output_dir
        self.progress_file = f"{progress_file}.csv"
        self.csv_file_path = csv_file_path
        self.completed_jobs = self._load_progress()

    def _load_progress(self):
        """
        Load completed jobs from progress file to avoid re-running.
        """
        if os.path.exists(self.progress_file):
            return set(pd.read_csv(self.progress_file)["job_id"])
        return set()

    def _save_progress(self, job_id: str):
        """
        Save the completed job to progress file.
        """
        with open(self.progress_file, "a") as f:
            f.write(f"{job_id}\n")

    def _generate_combinations(self):
        """
        Generate all single-monomer, two-monomer, and three-monomer combinations with solvents.
        """
        all_combinations = []
        for _, solvent in self.solvent_df.iterrows():
            solvent_name = solvent["name"]
            solvent_smiles = solvent["SMILES"]
            solvent_compressibility = float(solvent["compressibility"])
            solvent_density = float(solvent["density"])

            # Single-monomer cases
            for monomer in self.monomer_smiles:
                all_combinations.append(
                    (
                        [monomer],
                        solvent_name,
                        solvent_smiles,
                        solvent_density,
                        solvent_compressibility,
                    )
                )

            # Two-monomer copolymers (order-independent, no duplicates)
            for monomer_pair in itertools.combinations(self.monomer_smiles, 2):
                all_combinations.append(
                    (
                        list(monomer_pair),
                        solvent_name,
                        solvent_smiles,
                        solvent_density,
                        solvent_compressibility,
                    )
                )

            # Three-monomer copolymers
            for monomer_triplet in itertools.combinations(self.monomer_smiles, 3):
                all_combinations.append(
                    (
                        list(monomer_triplet),
                        solvent_name,
                        solvent_smiles,
                        solvent_density,
                        solvent_compressibility,
                    )
                )

        return all_combinations

    def run(self):
        """
        Iterates through all polymer-solvent-temperature combinations and executes the workflow.
        """
        combinations = self._generate_combinations()
        error_log = []

        for (
            monomer_list,
            solvent_name,
            solvent_smiles,
            solvent_density,
            solvent_compressibility,
        ) in combinations:
            for num_units in self.num_units:
                for temp in self.temperatures:
                    job_id = (
                        f"{'_'.join(monomer_list)}_{solvent_name}_{num_units}_{temp}"
                    )

                    if job_id in self.completed_jobs:
                        logger.info(f"Skipping {job_id}, already completed.")
                        continue

                    try:
                        logger.info(f"Running simulation for: {job_id}")

                        workflow = JoinedAtomisticPolymerWorkflow(
                            monomer_smiles=monomer_list,
                            num_units=num_units,
                            solvent_name=solvent_name,
                            solvent_smiles=solvent_smiles,
                            solvent_density=solvent_density,
                            temperatures=[temp],
                            output_dir=self.output_dir,
                            solvent_compressibility=solvent_compressibility,
                            csv_file_path=self.csv_file_path,
                        )

                        workflow.run()

                        # Log progress
                        self._save_progress(job_id)

                    except ValueError as e:
                        logger.error(f"ValueError in {job_id}: {e}")
                        error_log.append((job_id, str(e)))
                    except Exception as e:
                        logger.error(f"Unexpected error in {job_id}: {e}")
                        error_log.append((job_id, str(e)))

        # Save error log
        error_df = pd.DataFrame(error_log, columns=["Job ID", "Error"])
        error_df.to_csv("error_log.csv", index=False)
        logger.info("All simulations completed!")
