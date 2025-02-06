from modules.directory_parser.polymer_directory_parser import PolymerDirectoryParser
from modules.directory_parser.solvent_directory_parser import SolventDirectoryParser
from modules.workflows.separated.gromacs.joined import JoinedAtomisticPolymerWorkflow
from modules.utils.shared.file_utils import check_directory_exists
import logging
import os
import pandas as pd
from typing import List

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class SimulationManager:
    error_csv_filename = "error_log"
    output_csv_filename = "outputs"  # .csv

    def __init__(
        self,
        parameterised_poly_dir: str,
        parameterised_sol_dir: str,
        temperatures: List[float],
        output_dir: str,
        identifying_tag: str = "",
        solvent_start_idx: int = 0,
        solvent_end_idx: int = None,
        polymer_start_idx: int = 0,
        polymer_end_idx: int = None,
    ):
        self.parameterised_solvents = SolventDirectoryParser(
            parameterised_sol_dir
        ).parse_directory()
        self.parameterised_polymers = PolymerDirectoryParser(
            parameterised_poly_dir
        ).parse_directory()

        # truncating so that I can split HPC submissions into multiple jobs
        solvent_length = len(self.parameterised_solvents)
        if solvent_end_idx is None or solvent_end_idx > solvent_length:
            solvent_end_idx = solvent_length

        polymer_length = len(self.parameterised_polymers)
        if polymer_end_idx is None or polymer_end_idx > polymer_length:
            polymer_end_idx = polymer_length

        self.parameterised_solvents = self.parameterised_solvents[
            solvent_start_idx:solvent_end_idx
        ]
        self.parameterised_polymers = self.parameterised_polymers[
            polymer_start_idx:polymer_end_idx
        ]
        self.output_dir = output_dir
        check_directory_exists(output_dir, make_dirs=True)
        self.temperatures = temperatures
        self.error_file_path = os.path.join(
            output_dir, f"{self.error_csv_filename}_{identifying_tag}.csv"
        )

    def _get_sanitised_monomer_smiles(self, monomer_smiles: List[str]) -> str:
        return "_".join(monomer_smiles)

    def run(self):
        for (
            parameterised_solvent,
            solvent,
            solvent_smiles,
        ) in self.parameterised_solvents:
            for (
                parameterised_polymer,
                monomer_smiles,
                n_units,
            ) in self.parameterised_polymers:
                for temp in self.temperatures:

                    job_id = f"{solvent.name}_{self._get_sanitised_monomer_smiles(monomer_smiles=monomer_smiles)}_{n_units}_{temp}"
                    try:
                        workflow = JoinedAtomisticPolymerWorkflow(
                            parameterised_polymer=parameterised_polymer,
                            monomer_smiles=monomer_smiles,
                            num_units=n_units,
                            solvent=solvent,
                            solvent_smiles=solvent_smiles,
                            parameterised_solvent=parameterised_solvent,
                            temperature=temp,
                            output_dir=self.output_dir,
                            csv_file_path=self.output_csv_filename,
                        ).run()
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

        error_entry = pd.DataFrame(
            [[job_id, error_message]], columns=["Job ID", "Error"]
        )
        error_entry.to_csv(self.error_file_path, mode="a", header=False, index=False)
