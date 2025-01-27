from open_mscg_topol_gen import OpenMSCGTopolGenerator
from itertools import combinations, product
from typing import List, Dict, Tuple, Optional, Set, Union
from mscg.cli import cgfm
from modules.shared.utils.file_utils import check_directory_exists, check_file_type
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class OpenMSCGForceMatcher:
    non_bonded_interaction_params = {
        "poly:poly": {"min": 0.0, "max": 15.0, "res": 0.5},
        "sol:poly": {"min": 0.0, "max": 15.0, "res": 0.3},
        "sol:sol": {"min": 0.0, "max": 10.0, "res": 0.2},
        "ion:ion": {"min": 0.0, "max": 8.0, "res": 0.1},
        "ion:sol": {"min": 0.0, "max": 10.0, "res": 0.2},
        "ion:poly": {"min": 0.0, "max": 12.0, "res": 0.2},
    }
    bonded_interaction_params = {"min": 1.0, "max": 2.5, "res": 0.01}
    cut = 15.0

    def __init__(
        self,
        topol_generator: OpenMSCGTopolGenerator,
        traj_path: str,
        sol_name: str = "SOL",
        ions_list: List[str] = ["NA", "CL"],
        model="BSpline",
    ):

        self._model = model
        self._check_input_file_types(
            top_path=topol_generator.top_path, traj_path=traj_path
        )
        self._traj_path = traj_path
        self._top_path = topol_generator.top_path
        self.mols_list = topol_generator.mols_list
        self.ions_list = ions_list
        self.sol_name = sol_name

        self.all_non_bonded_pairs = None
        self.all_bonded_pairs = None
        self._validate_params(
            cut=self.cut,
            non_bonded_interaction_params=self.non_bonded_interaction_params,
            bonded_interaction_params=self.bonded_interaction_params,
        )

        self._generate_all_pairs()

    def _check_input_file_types(
        self, top_path: Optional[str] = None, traj_path: Optional[str] = None
    ):
        if top_path:
            check_file_type(file_path=top_path, expected_file_type="top")
        if traj_path:
            check_file_type(file_path=traj_path, expected_file_type="lammpstrj")

    def _validate_params(
        self,
        cut: float,
        non_bonded_interaction_params: Dict[str, Dict[str, float]],
        bonded_interaction_params: Dict[str, float],
    ):
        max_value = max(
            max(entry["max"] for entry in non_bonded_interaction_params.values()),
            bonded_interaction_params["max"],
        )
        if cut < max_value:
            raise ValueError(
                f"cut value {cut} must be >= max value in interaction params {max_value} (as cut value is the absolute max cut-off val)"
            )

    @property
    def model(self) -> str:
        return self._model

    @model.setter
    def model(self, model: str) -> None:
        self._model = model

    @property
    def cut(self) -> float:
        return self._cut

    @cut.setter
    def cut(self, cut: float) -> None:
        self._validate_params(
            cut=cut,
            bonded_interaction_params=self.bonded_interaction_params,
            non_bonded_interaction_params=self.non_bonded_interaction_params,
        )
        self._cut = cut

    @property
    def traj_path(self) -> str:
        return self._traj_path

    @traj_path.setter
    def traj_path(self, traj_path: str) -> None:
        self._check_input_file_types(traj_path=traj_path)
        self._traj_path = traj_path

    @property
    def top_path(self) -> str:
        return self._top_path

    @top_path.setter
    def top_path(self, top_path: str) -> None:
        self._check_input_file_types(top_path=top_path)
        self._top_path = top_path

    def _is_multiple(self):
        for mol in self.mols_list:
            if mol["repeat_count"] > 1:
                mol["is_multiple"] = True
            else:
                mol["is_multiple"] = False

    def configure_non_bonded_parameters(
        self,
        interaction_type: Optional[str],
        min: Optional[float] = None,
        max: Optional[float] = None,
        resolution: Optional[float] = None,
    ) -> None:
        """
        Configures paramaters for the interaction type. Inputting None as an entry will leave it as the default params. Default params are the following:
            "poly:poly": (min: 0.0, max: 15.0, res: 0.5)
            "sol:poly":  (min: 0.0, max: 15.0, res: 0.3)
            "sol:sol":   (min: 0.0, max: 10.0, res: 0.2)
            "ion:ion":   (min: 0.0, max: 8.0, res: 0.1)
            "ion:sol":   (min: 0.0, max: 10.0, res: 0.2)
            "ion:poly":  (min: 0.0, max: 12.0, res: 0.2)

        :param interaction_type: interaction type to configure settings for, set to None to configure all. Allowed inputs are: "poly:poly, sol:poly, sol:sol, ion:ion, ion:sol, ion:poly."
        :type interaction_type: str
        :param min: The minimum separation distance for the interaction to be calculated for in Angstroms, defaults to None
        :type min: Optional[str], optional
        :param max: The maximum separation distance for the interaction to be calculated for in Angstroms, defaults to None
        :type max: Optional[str], optional
        :param resolution: the size of the bin spacing for BSpline interactions, defaults to None
        :type resolution: Optional[str], optional
        :raises ValueError: _description_
        :raises ValueError: _description_
        :return: None (sets self.interaction_params)
        :rtype: None
        """
        valid_interactions = set(self.non_bonded_interaction_params.keys())

        if interaction_type and interaction_type not in valid_interactions:
            raise ValueError(
                f"Invalid interaction type '{interaction_type}'. Allowed values are: {', '.join(valid_interactions)}"
            )

        if max:
            self._validate_params(
                cut=self.cut,
                non_bonded_interaction_params={"max": {"max": max}},
                bonded_interaction_params=self.bonded_interaction_params,
            )

        interactions_to_update = (
            [interaction_type]
            if interaction_type
            else list(self.non_bonded_interaction_params.keys())
        )

        for interaction in interactions_to_update:
            current_params = self.non_bonded_interaction_params[interaction]
            self.non_bonded_interaction_params[interaction] = {
                "min": min if min is not None else current_params["min"],
                "max": max if max is not None else current_params["max"],
                "res": resolution if resolution is not None else current_params["res"],
            }

        logging.info(
            f"Updated interaction parameters: {self.non_bonded_interaction_params}"
        )

    def configure_bonded_parameters(
        self,
        min: Optional[float] = None,
        max: Optional[float] = None,
        resolution: Optional[float] = None,
    ) -> None:
        """
        Configures paramaters for the interaction type for bonded molecule.
        Default is the following:
        (min: 1.0, max: 2.5, res: 0.0)
        :param interaction_type: interaction type to configure settings for, set to None to configure all. Allowed inputs are: "poly:poly, sol:poly, sol:sol, ion:ion, ion:sol, ion:poly."
        :type interaction_type: str
        :param min: The minimum separation distance for the interaction to be calculated for in Angstroms, defaults to None
        :type min: Optional[str], optional
        :param max: The maximum separation distance for the interaction to be calculated for in Angstroms, defaults to None
        :type max: Optional[str], optional
        :param resolution: the size of the bin spacing for BSpline interactions, defaults to None
        :type resolution: Optional[str], optional
        :raises ValueError: _description_
        :raises ValueError: _description_
        :return: None (sets self.interaction_params)
        :rtype: None
        """

        current_params = self.bonded_interaction_params
        if max:
            self._validate_params(
                cut=self.cut,
                non_bonded_interaction_params=self.non_bonded_interaction_params,
                bonded_interaction_params={"max": max},
            )
        self.bonded_interaction_params = {
            "min": min if min is not None else current_params["min"],
            "max": max if max is not None else current_params["max"],
            "res": resolution if resolution is not None else current_params["res"],
        }

        logging.info(
            f"Updated interaction parameters: {self.non_bonded_interaction_params}"
        )

    def _determine_mol_type(self, site: str) -> str:
        if site in self.ions_list:
            return "ion"
        elif site == self.sol_name:
            return "sol"
        else:
            return "poly"

    def _determine_interation_type(self, sites_pair: str) -> str:
        """
        Gets the interaction type for a given site pair.

        :param site_pair: sites pair to get interaction for, in format: "site1:site2"
        :type site_pair: str
        :return: site type matching self.interaction_params key definitions
        :rtype: str
        """
        site1, site2 = sites_pair.split(":")
        types = {self._determine_mol_type(site1), self._determine_mol_type(site2)}

        if "ion" in types:
            if "sol" in types:
                return "ion:sol"
            if "poly" in types:
                return "ion:poly"
            return "ion:ion"
        if "sol" in types:
            return "sol:poly" if "poly" in types else "sol:sol"
        return "poly:poly"

    def _get_molecule_non_bonded_pairs(
        self, molecule: Dict[str, Union[int, List[List[Union[str, int]]]]]
    ) -> List[str]:
        include_self = molecule["is_multiple"]
        cg_sites = [cg_site[0] for cg_site in molecule["sites"]]
        if include_self:
            site_pairs = set(product(cg_sites, repeat=2))
        else:
            site_pairs = set(combinations(cg_sites, 2))

        unique_pairs = set(tuple(sorted(site_pair)) for site_pair in site_pairs)

        return list(unique_pairs)

    def _get_molecule_bonded_pairs(
        self, molecule: Dict[str, Union[int, List[List[Union[str, int]]]]]
    ) -> List[str]:
        cg_bonded_sites = set()
        molecule_sites = molecule["sites"]
        if len(molecule_sites) < 2:
            return cg_bonded_sites

        for i in range(len(molecule_sites) - 1):
            bond = (molecule_sites[i][0], molecule_sites[i + 1][0])
            cg_bonded_sites.add(tuple(sorted(bond)))
        return cg_bonded_sites

    def _format_pairs_set(self, cg_pairs_set: Set[str]) -> List[str]:
        return [":".join(sorted(cg_pair)) for cg_pair in cg_pairs_set]

    def _generate_all_pairs(self) -> None:
        self._is_multiple()
        all_non_bonded_pairs = set()
        all_bonded_pairs = set()

        for molecule in self.mols_list:
            non_bonded_pairs = self._get_molecule_non_bonded_pairs(molecule)
            all_non_bonded_pairs.update(non_bonded_pairs)

            bonded_pairs = self._get_molecule_bonded_pairs(molecule)
            all_bonded_pairs.update(bonded_pairs)

        self.all_non_bonded_pairs = self._format_pairs_set(all_non_bonded_pairs)
        self.all_bonded_pairs = self._format_pairs_set(all_bonded_pairs)

    def _get_non_bonded_params(
        self, sites_pair: str, include_ions: bool = False
    ) -> str:
        interaction_type = self._determine_interation_type(sites_pair)
        if not include_ions and "ion" in interaction_type:
            return None

        params = self.non_bonded_interaction_params[interaction_type]
        return params

    def _format_pair_input(self, params: Dict[str, float], sites_pair: str) -> str:
        return f"model={self.model},type={sites_pair},min={params['min']},max={params['max']},resolution={params['res']}"

    def _create_non_bonded_pair_inputs(
        self, include_ions: bool = False
    ) -> Tuple[List[str], List[str]]:
        logger.info(
            f"Generating non bonded pairs inputs with parameters: {self.non_bonded_interaction_params}"
        )
        non_bonded_inputs = []
        for non_bonded_pair in self.all_non_bonded_pairs:
            non_bonded_input = self._get_non_bonded_params(
                non_bonded_pair, include_ions
            )
            if non_bonded_input:
                non_bonded_inputs.append(
                    self._format_pair_input(
                        params=non_bonded_input, sites_pair=non_bonded_pair
                    )
                )
        return non_bonded_inputs

    def _create_bonded_pair_inputs(self) -> Tuple[List[str], List[str]]:
        logger.info(
            f"Generating bonded pairs inputs with parameters: {self.bonded_interaction_params}"
        )
        bonded_params = self.bonded_interaction_params
        bonded_inputs = []
        for bonded_pair in self.all_bonded_pairs:
            if bonded_pair:
                bonded_inputs.append(
                    self._format_pair_input(
                        params=bonded_params, sites_pair=bonded_pair
                    )
                )
        return bonded_inputs

    def run(
        self,
        file_name: str,
        output_dir: Optional[str] = None,
        include_ions: bool = False,
    ):

        non_bonded_input = self._create_non_bonded_pair_inputs(include_ions)
        bonded_inputs = self._create_bonded_pair_inputs()
        logger.info(
            f"Using cgfm with cutoff: {self.cut}, non-bonded inputs: {non_bonded_input}, bonded inputs: {bonded_inputs}"
        )
        if output_dir:
            check_directory_exists(output_dir)
            output_path = f"{output_dir}/{file_name}"
        else:
            output_path = f"{file_name}"

        logger.info(
            f"Running cgfm with inputs: top: {self.top_path}, traj: {self.traj_path}"
        )
        logger.info(f"Output will be saved to: {output_path}.g")

        cgfm.main(
            top=self.top_path,
            traj=self.traj_path,
            cut=self.cut,
            save=output_path,
            pair=non_bonded_input,
            bond=bonded_inputs,
        )
