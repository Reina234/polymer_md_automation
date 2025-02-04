import subprocess
import numpy as np
import csv
import os

########################################################################
# 1. OUTPUT CONTAINER
########################################################################


class GromacsOutputs:
    """
    Holds file paths to GROMACS output files: .tpr, .xtc, .gro, etc.
    """

    def __init__(self, tpr="simulation.tpr", xtc="trajectory.xtc", gro="final.gro"):
        self.tpr = tpr
        self.xtc = xtc
        self.gro = gro


########################################################################
# 2. INDEX MANAGER
########################################################################


class GromacsAnalysis:
    """
    Performs various GROMACS analyses using a provided GromacsOutputs
    object and an index file. Ensures we select the correct group (e.g., "Polymer").
    """

    def __init__(self, outputs: GromacsOutputs, index_file="index.ndx"):
        self.outputs = outputs
        self.index_file = index_file

    def extract_radius_of_gyration(self):
        """
        Uses 'gmx gyrate' to compute radius of gyration for the "Polymer" group.
        """
        cmd = [
            "gmx",
            "gyrate",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            "gyrate.xvg",
        ]
        # We must pick the "Polymer" group from the index file:
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate(
            "Polymer\n"
        )  # The exact group name must match what we used in make_ndx
        p.wait()

        data = np.loadtxt("gyrate.xvg", comments=("#", "@"))
        Rg_mean = np.mean(data[:, 1])
        Rg_std = np.std(data[:, 1])
        return Rg_mean, Rg_std

    def extract_diffusion_coefficient(self):
        """
        Uses 'gmx msd' to compute mean-squared displacement for the "Polymer" group.
        """
        cmd = [
            "gmx",
            "msd",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            "msd.xvg",
        ]
        # Select the polymer group
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt("msd.xvg", comments=("#", "@"))
        t = data[:, 0]
        msd = data[:, 1]
        start = int(0.8 * len(t))
        slope, _ = np.polyfit(t[start:], msd[start:], 1)
        D = slope / 6.0  # 3D diffusion
        return D

    def extract_sasa(self):
        """
        Uses 'gmx sasa' to compute SASA for the "Polymer" group.
        """
        cmd = [
            "gmx",
            "sasa",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            "sasa.xvg",
        ]
        # Select the polymer group
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt("sasa.xvg", comments=("#", "@"))
        SASA_mean = np.mean(data[:, 1])
        SASA_std = np.std(data[:, 1])
        return SASA_mean, SASA_std

    def extract_end_to_end_distance(self):
        """
        Uses 'gmx distance' to compute end-to-end distance for the "Polymer" group.
        Typically you'd specify two groups or atoms, but here's a placeholder approach.
        """
        cmd = [
            "gmx",
            "distance",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-o",
            "end_to_end.xvg",
        ]
        # This might require specifying pairs of atoms. If 'distance' prompts,
        # you'll feed in e.g. "Polymer" or specific selection.
        # In a real system, you might define special groups for "first_atom" vs. "last_atom".
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        # Example: If we had created groups for "Polymer_end1" and "Polymer_end2"
        # p.communicate("Polymer_end1\nPolymer_end2\n")
        # For now, just a placeholder:
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt("end_to_end.xvg", comments=("#", "@"))
        E2E_mean = np.mean(data[:, 1])
        E2E_std = np.std(data[:, 1])
        return E2E_mean, E2E_std

    # Add any additional extractions (shape descriptors, RDF, persistence length, etc.)
    # using similar patterns with -n index.ndx and group selections.


########################################################################
# 4. MAIN WORKFLOW
########################################################################


class SimulationWorkflow:
    """
    Orchestrates:
      - Running or calling the simulation (dummy here)
      - Creating the index file
      - Extracting desired mesoscale properties
      - Saving to CSV
    """

    def __init__(
        self,
        solvent_smiles,
        compressibility,
        monomer_smiles_list,
        N,
        T,
        csv_filename="results.csv",
    ):
        self.solvent_smiles = solvent_smiles
        self.compressibility = compressibility
        self.monomer_smiles_list = monomer_smiles_list
        self.N = N
        self.T = T
        self.csv_filename = csv_filename

    def run_simulation(self):
        """
        Dummy function that 'runs' the simulation and returns a GromacsOutputs object.
        Replace with actual GROMACS calls or submission scripts.
        """
        return GromacsOutputs(
            tpr="simulation.tpr", xtc="trajectory.xtc", gro="final.gro"
        )

    def run_workflow(self):
        """
        Full workflow:
          1) run simulation,
          2) create index file,
          3) extract properties,
          4) write CSV row
        """
        # 1) Run the simulation
        outputs = self.run_simulation()

        # 2) Create index file for analysis
        index_manager = GromacsIndexManager(outputs.gro)
        index_manager.create_index_file()

        # 3) Perform analyses
        analyzer = GromacsAnalysis(outputs, index_file=index_manager.index_file)
        Rg_mean, Rg_std = analyzer.extract_radius_of_gyration()
        D = analyzer.extract_diffusion_coefficient()
        SASA_mean, SASA_std = analyzer.extract_sasa()
        E2E_mean, E2E_std = analyzer.extract_end_to_end_distance()
        # etc. (Add more as desired...)

        # 4) Prepare a dictionary for CSV output
        row_data = {
            "solvent_smiles": self.solvent_smiles,
            "compressibility": self.compressibility,
            "monomer_smiles_list": ";".join(self.monomer_smiles_list),
            "N": self.N,
            "T": self.T,
            "Rg_mean": Rg_mean,
            "Rg_std": Rg_std,
            "Diffusion_Coefficient": D,
            "SASA_mean": SASA_mean,
            "SASA_std": SASA_std,
            "E2E_mean": E2E_mean,
            "E2E_std": E2E_std,
        }

        # 5) Append to CSV
        file_exists = os.path.isfile(self.csv_filename)
        with open(self.csv_filename, mode="a", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=row_data.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(row_data)

        print(f"Row appended to {self.csv_filename}")

        return row_data  # Optional: return data if needed


########################################################################
# EXAMPLE USAGE
########################################################################

if __name__ == "__main__":
    # Example input
    solvent_smiles = "O"
    compressibility = 4.5e-10
    monomer_smiles_list = ["C=C"]
    N = 100
    T = 300
    csv_filename = "sim_results.csv"

    # Run the workflow
    workflow = SimulationWorkflow(
        solvent_smiles, compressibility, monomer_smiles_list, N, T, csv_filename
    )
    workflow.run_workflow()
