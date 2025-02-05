import itertools
import logging
import os
import pandas as pd
from typing import List, Tuple
from modules.workflows.atomistic.no_parametiser_joined import (
    NoJoinedAtomisticPolymerWorkflow,
)
from modules.workflows.separated.parametiser.polymer import PolymerParametiser
from modules.workflows.separated.parametiser.solvent import SolventParametiser
from modules.utils.shared.file_utils import check_directory_exists

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PolymerPreparationManager:
    def __init__(
        self,
        monomer_smiles: List[str],
        output_dir: str,
        num_units: List[int] = [5, 10, 20],
    ):
        """
        Initializes the simulation manager.
        :param solvent_csv: Path to the CSV file containing solvent properties.
        :param monomer_smiles: List of monomer SMILES strings.
        :param output_dir: Directory where results should be stored.
        :param progress_file: File for tracking progress.
        """
        self.monomer_smiles = monomer_smiles
        self.num_units = num_units
        self.output_dir = output_dir

    def _generate_combinations(self):
        """
        Generate all single-monomer, two-monomer, and three-monomer combinations with solvents.
        """
        all_combinations = []
        for monomer in self.monomer_smiles:
            all_combinations.append([monomer])

        # Two-monomer copolymers (order-independent, no duplicates)
        for monomer_pair in itertools.combinations(self.monomer_smiles, 2):
            all_combinations.append(list(monomer_pair))

        # Three-monomer copolymers
        for monomer_triplet in itertools.combinations(self.monomer_smiles, 3):
            all_combinations.append(list(monomer_triplet))

        return all_combinations

    def run(self):
        monomer_combinations = self._generate_combinations()
        for monomer_smiles in monomer_combinations:
            for num_units in self.num_units:
                generator = PolymerParametiser(
                    monomer_smiles, num_units, output_dir=self.output_dir
                )
                generator.run()


class SolventPreparationManager:
    def __init__(
        self,
        solvent_csv: str,
        output_dir: str,
        verbose: bool = False,
    ):
        """
        Initializes the simulation manager.
        :param solvent_csv: Path to the CSV file containing solvent properties.
        :param monomer_smiles: List of monomer SMILES strings.
        :param output_dir: Directory where results should be stored.
        :param progress_file: File for tracking progress.
        """
        check_directory_exists(output_dir, make_dirs=True)
        self.verbose = verbose
        self.output_dir = output_dir
        self.solvent_df = pd.read_csv(solvent_csv)

    def run(self):

        for _, solvent in self.solvent_df.iterrows():
            solvent_name = solvent["name"]
            solvent_smiles = solvent["SMILES"]
            solvent_compressibility = float(solvent["compressibility"])
            solvent_density = float(solvent["density"])
            solvent_parametiser = SolventParametiser(
                solvent_name=solvent_name,
                solvent_smiles=solvent_smiles,
                solvent_compressibility=solvent_compressibility,
                solvent_density=solvent_density,
                output_dir=self.output_dir,
                verbose=self.verbose,
            )
            solvent_parametiser.run()


class GromacsManager:
    def __init__(
        self,
        solvent_csv: str,
        input_dir: str,
        output_dir: str,
        csv_file: str = "outputs",
        temperatures: List[int] = [280, 298, 346],
        error_log: str = "error_log",
    ):
        """
        Initializes the simulation manager.
        :param solvent_csv: Path to the CSV file containing solvent properties.
        :param monomer_smiles: List of monomer SMILES strings.
        :param output_dir: Directory where results should be stored.
        :param progress_file: File for tracking progress.
        """
        check_directory_exists(output_dir, make_dirs=True)
        check_directory_exists(input_dir, make_dirs=False)
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.solvent_df = pd.read_csv(solvent_csv)
        self.temperatures = temperatures
        self.csv_file_path = os.path.join(output_dir, f"{csv_file}.csv")
        self.error_log_file = os.path.join(output_dir, f"{error_log}.csv")

    def run(self):
        error_log = []

        for _, solvent in self.solvent_df.iterrows():
            solvent_name = solvent["name"]
            solvent_smiles = solvent["SMILES"]
            solvent_compressibility = float(solvent["compressibility"])
            solvent_density = float(solvent["density"])
            job_id = f"{solvent_name}_{solvent_smiles}"
            for temp in self.temperatures:
                try:
                    workflow = NoJoinedAtomisticPolymerWorkflow(
                        input_dir=self.input_dir,
                        solvent_name=solvent_name,
                        solvent_smiles=solvent_smiles,
                        solvent_density=solvent_density,
                        temperatures=[temp],
                        output_dir=self.output_dir,
                        solvent_compressibility=solvent_compressibility,
                        csv_file_path=self.csv_file_path,
                    )

                    workflow.run()
                except ValueError as e:
                    error_message = f"ValueError in {job_id}: {e}"
                    logger.error(error_message)
                    self._log_error(job_id, str(e))

                except Exception as e:
                    error_message = f"Unexpected error in {job_id}: {e}"
                    logger.error(error_message)
                    self._log_error(job_id, str(e))

        logger.info("All simulations completed!")

    def _log_error(self, job_id: str, error_message: str):
        """
        Logs an error immediately to the error log file.

        :param job_id: Identifier for the job.
        :param error_message: Description of the error.
        """
        error_entry = pd.DataFrame(
            [[job_id, error_message]], columns=["Job ID", "Error"]
        )
        error_entry.to_csv(self.error_log_file, mode="a", header=False, index=False)
