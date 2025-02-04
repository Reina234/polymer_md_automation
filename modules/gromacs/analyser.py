from config.data_models.output_types import GromacsOutputs
from modules.gromacs.index_manager import GromacsIndexManager
from config.paths import TEMP_DIR
import subprocess
from typing import List
import os
import numpy as np


class GromacsAnalyser:
    RADIUS_GYRATION_FILE = "gyrate.xvg"
    DIFFUSION_COEFFICIENT_FILE = "msd.xvg"
    SASA_FILE = "sasa.xvg"
    END_TO_END_DISTANCE_FILE = "end_to_end.xvg"

    def __init__(
        self,
        outputs: GromacsOutputs,
        poly_resname: str = "UNL",
        ion_resnames: List[str] = ["NA", "CL"],
        output_dir: str = TEMP_DIR,
    ):
        self.outputs = outputs
        print("!!!!!!!!!!1")
        print(self.outputs)
        self.output_dir = output_dir
        os.makedirs(self.output_dir, exist_ok=True)

        self.index_file = self._create_index_file(
            gro_file=outputs.gro,
            output_dir=self.output_dir,
            poly_resname=poly_resname,
            ion_resnames=ion_resnames,
        )

    def _create_index_file(
        self, gro_file: str, output_dir: str, poly_resname: str, ion_resnames: List[str]
    ) -> str:
        index_manager = GromacsIndexManager(
            gro_file=gro_file,
            output_dir=output_dir,
            poly_name=poly_resname,
            ion_names=ion_resnames,
        )
        return index_manager.create_index_file()

    def _get_output_path(self, basename: str) -> str:
        """
        Constructs the full output path inside output_dir.
        """
        return os.path.join(self.output_dir, basename)

    def extract_radius_of_gyration(self):
        output_path = self._get_output_path(self.RADIUS_GYRATION_FILE)

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
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt(output_path, comments=("#", "@"))
        Rg_mean = np.mean(data[:, 1])
        Rg_std = np.std(data[:, 1])
        return Rg_mean, Rg_std

    def extract_diffusion_coefficient(self):
        output_path = self._get_output_path(self.DIFFUSION_COEFFICIENT_FILE)

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
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt(output_path, comments=("#", "@"))
        t = data[:, 0]
        msd = data[:, 1]
        start = int(0.8 * len(t))
        slope, _ = np.polyfit(t[start:], msd[start:], 1)
        D = slope / 6.0
        return D

    def extract_sasa(self):
        output_path = self._get_output_path(self.SASA_FILE)

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
            output_path,
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)
        p.communicate("Polymer\n")
        p.wait()

        data = np.loadtxt(output_path, comments=("#", "@"))
        SASA_mean = np.mean(data[:, 1])
        SASA_std = np.std(data[:, 1])
        return SASA_mean, SASA_std

    def extract_end_to_end_distance(self):
        """
        Computes the end-to-end distance of the polymer using 'gmx distance'.
        Uses the first and last polymer carbon atoms as reference points.
        """
        output_path = self._get_output_path(self.END_TO_END_DISTANCE_FILE)

        cmd = [
            "gmx",
            "distance",
            "-f",
            self.outputs.xtc,
            "-s",
            self.outputs.tpr,
            "-n",
            self.index_file,
            "-oall",
            output_path,  # FIX: use -oall instead of -o
        ]

        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, text=True)

        # FIX: Ensure correct group selection
        selection_input = f"Polymer_Carbons_Start_and_End\n"
        p.communicate(input=selection_input)
        p.wait()

        if not os.path.exists(output_path):
            raise FileNotFoundError(f"Output file not found: {output_path}")

        # Load and process the output file
        try:
            data = np.loadtxt(output_path, comments=("#", "@"))
            E2E_mean = np.mean(data[:, 1])
            E2E_std = np.std(data[:, 1])
            return E2E_mean, E2E_std
        except Exception as e:
            raise RuntimeError(f"Error reading end-to-end distance file: {e}")
