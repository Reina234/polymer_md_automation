import os
import subprocess
from preprocessing.metadata_tracker import MetadataTracker
from config.paths import (
    MDP_FULL_PATHS,
    TemplatedMdps,
    BASE_OUTPUT_DIR,
    GROMACS_OUTPUT_SUBDIR,
)
from gromacs.base_gromacs_command import BaseGromacsCommand
from typing import Optional, List, Tuple


class IonAdder(BaseGromacsCommand):
    OUTPUT_NAME = "polymer_neutralized.gro"
    OUTPUT_TPR_NAME = "ions.tpr"

    def __init__(self, metadata_tracker: MetadataTracker):
        self.metadata_tracker = metadata_tracker

    def run(
        self,
        solvated_gro_path: str,
        input_topol_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        minim_mdp_path: Optional[str] = None,
    ):
        grommp_command, output_tpr_path, output_dir = self._create_grommp_command(
            solvated_gro_path=solvated_gro_path,
            input_topol_path=input_topol_path,
            run_name=run_name,
            output_base_dir=output_base_dir,
            minim_mdp_path=minim_mdp_path,
        )
        self._execute(grommp_command)

        additional_command, output_gro_path = self._create_genion_command(
            output_tpr_path=output_tpr_path,
            output_dir=output_dir,
            input_topol_path=input_topol_path,
        )
        self._execute(additional_command)
        if self.metadata_tracker:
            self._update_metadata(
                solvated_gro_path=solvated_gro_path,
                input_topol_path=input_topol_path,
                run_name=run_name,
                output_base_dir=output_base_dir,
                minim_mdp_path=minim_mdp_path,
            )

        return output_gro_path

    def _create_grommp_command(
        self,
        solvated_gro_path: str,
        input_topol_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        minim_mdp_path: Optional[str] = None,
        additional_notes: Optional[str] = None,
    ) -> Tuple[List, str, str]:
        output_dir = os.path.join(output_base_dir, run_name, GROMACS_OUTPUT_SUBDIR)
        os.makedirs(output_dir, exist_ok=True)

        if not minim_mdp_path:
            minim_mdp_path = MDP_FULL_PATHS[TemplatedMdps.MINIM.value]

        output_tpr_path = os.path.join(output_dir, self.OUTPUT_TPR_NAME)

        grompp_command = [
            "gmx",
            "grompp",
            "-f",
            minim_mdp_path,
            "-c",
            solvated_gro_path,
            "-p",
            input_topol_path,
            "-o",
            output_tpr_path,
            "-maxwarn",
            "1",
        ]
        return grompp_command, output_tpr_path, output_dir

    def _create_genion_command(
        self,
        output_tpr_path: str,
        output_dir: str,
        input_topol_path: str,
    ) -> Tuple[List, str]:
        output_gro_path = os.path.join(output_dir, self.OUTPUT_NAME)
        genion_command = [
            "gmx",
            "genion",
            "-s",
            output_tpr_path,
            "-o",
            output_gro_path,
            "-p",
            input_topol_path,
            "-pname",
            "NA",
            "-nname",
            "CL",
            "-neutral",
        ]

        return genion_command, output_gro_path

    def metadata(
        self,
        solvated_gro_path: str,
        input_topol_path: str,
        run_name: str,
        output_base_dir: str = BASE_OUTPUT_DIR,
        minim_mdp_path: Optional[str] = None,
        additional_notes: Optional[str] = None,
    ):
        return {
            "program(s) used": "GROMACS solvate",
            "details": f"added ions to {solvated_gro_path} to neutralised",
            "action(s)": f"saved at {output_base_dir}/{run_name}/{GROMACS_OUTPUT_SUBDIR}/{self.OUTPUT_NAME}",
            "additional_notes": additional_notes,
        }
