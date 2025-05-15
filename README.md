# Automated Polymer Simulations

This README is a work in progress.

This repo contains code for automating GROMACS molecular dynamics simulation, currently capable of handling addition homopolymers and simple alternating copolymers. The repository also has methods for automatically deriving coarse-grained force fields via force-matching methodology against short atomistic trajectories, as well as tools for constructing the full coarse-grained system. However, this has not yet been integrated into the full simulation automation due to time constraints, and exists just as manual-use functionality.

Overall workflow handling is done via the simulation manager, it uses CSV data for solvents information and a list for the monomer smiles (generates all combinations of 1, 2, 3 monomer type composition). This is called in main.py. Caching, error logging, output formatting is all handled by `SimulationManager`. Cache, logs, and temp folders will be automatically generated.

`parametiser_workflow.py` and `directory_parser` hold defunct methods that does the whole parameterisation workflow. This was part of initial separation attempts of the parameterisation workflow vs MD/GROMACS workflow for HPC uploading purposes. Example outputs, progress, and error CSVs are given in the samples folder.

All path imports assume the repo folder is top-level.

## General structure

```
📦repo
 ┣ 📂config
 ┃ ┣ 📂data_models
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜output_types.py
 ┃ ┃ ┗ 📜solvent.py
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜acpype_config.py
 ┃ ┣ 📜constants.py
 ┃ ┣ 📜mdp_workflow_config.py
 ┃ ┗ 📜paths.py
 ┣ 📂input_data
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜monomer_smiles.py
 ┃ ┣ 📜solvent_data.csv
 ┃ ┗ 📜solvent_data.xlsx
 ┣ 📂modules
 ┃ ┣ 📂acpype
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜acpype_parametizer.py
 ┃ ┃ ┗ 📜acpype_utils.py
 ┃ ┣ 📂cache_store
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_cache.py
 ┃ ┃ ┣ 📜equilibriated_atomistic_polymer_cache.py
 ┃ ┃ ┣ 📜file_cache.py
 ┃ ┃ ┣ 📜mdp_cache.py
 ┃ ┃ ┣ 📜pickle_cache.py
 ┃ ┃ ┗ 📜solvent_cache.py
 ┃ ┣ 📂cg_mappers
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_map_generator.py
 ┃ ┃ ┣ 📜martini_index_generator.py
 ┃ ┃ ┣ 📜martini_map_generator.py
 ┃ ┃ ┣ 📜multimol_map_generator.py
 ┃ ┃ ┣ 📜open_mscg_map_generator.py
 ┃ ┃ ┣ 📜pycgtool_map_generator.py
 ┃ ┃ ┗ 📜votca_map_generator.py
 ┃ ┣ 📂file_conversion
 ┃ ┃ ┣ 📂converters
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜base_converter.py
 ┃ ┃ ┃ ┣ 📜editconf_gro_to_pdb.py
 ┃ ┃ ┃ ┣ 📜editconf_pdb_to_gro.py
 ┃ ┃ ┃ ┗ 📜obabel_pdb_to_mol2_converter.py
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┗ 📜converter_factory.py
 ┃ ┣ 📂gromacs
 ┃ ┃ ┣ 📂equilibriation
 ┃ ┃ ┃ ┣ 📜base_workflow_step.py
 ┃ ┃ ┃ ┗ 📜full_equilibriation_workflow.py
 ┃ ┃ ┣ 📂parsers
 ┃ ┃ ┃ ┣ 📂data_models
 ┃ ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┃ ┗ 📜section.py
 ┃ ┃ ┃ ┣ 📂handlers
 ┃ ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┃ ┣ 📜data_handler.py
 ┃ ┃ ┃ ┃ ┗ 📜gro_handler.py
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜gromacs_parser.py
 ┃ ┃ ┃ ┗ 📜itp_parser.py
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜analyser.py
 ┃ ┃ ┗ 📜index_manager.py
 ┃ ┣ 📂lammps
 ┃ ┃ ┣ 📂parsers
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┗ 📜open_mscg_data_parser.py
 ┃ ┃ ┗ 📜__init__.py
 ┃ ┣ 📂moltemplate
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_molecule.py
 ┃ ┃ ┣ 📜moltemplate_system.py
 ┃ ┃ ┣ 📜moltemplate_utils.py
 ┃ ┃ ┣ 📜multimol_solvent.py
 ┃ ┃ ┣ 📜polymer.py
 ┃ ┃ ┗ 📜solvent.py
 ┃ ┣ 📂open_mscg
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜force_matcher.py
 ┃ ┃ ┣ 📜multimol_topol_exporter.py
 ┃ ┃ ┣ 📜multimol_traj_mapper.py
 ┃ ┃ ┣ 📜topol_exporter.py
 ┃ ┃ ┣ 📜topol_generator.py
 ┃ ┃ ┗ 📜trajectory_mapper.py
 ┃ ┣ 📂packmol
 ┃ ┃ ┣ 📂templates
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┗ 📜solvent_box_template.inp
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_packmol_operation.py
 ┃ ┃ ┗ 📜solvent_box.py
 ┃ ┣ 📂rdkit
 ┃ ┃ ┣ 📂polymer_builders
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜alternating_copolymer.py
 ┃ ┃ ┃ ┣ 📜base_polymer_generator.py
 ┃ ┃ ┃ ┗ 📜homopolymer_generator.py
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_molecule_generator.py
 ┃ ┃ ┣ 📜polymer_itp_scaler.py
 ┃ ┃ ┗ 📜solvent_generator.py
 ┃ ┣ 📂utils
 ┃ ┃ ┣ 📂atomistic
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜file_utils.py
 ┃ ┃ ┃ ┗ 📜mdp_utils.py
 ┃ ┃ ┣ 📂shared
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜calculation_utils.py
 ┃ ┃ ┃ ┣ 📜dataframe_utils.py
 ┃ ┃ ┃ ┗ 📜file_utils.py
 ┃ ┃ ┗ 📜__init__.py
 ┃ ┣ 📂workflows
 ┃ ┃ ┣ 📂atomistic
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜joined_workflow.py
 ┃ ┃ ┃ ┣ 📜polymer_equilibriator.py
 ┃ ┃ ┃ ┣ 📜polymer_parametizer.py
 ┃ ┃ ┃ ┗ 📜solvent_equilibriator.py
 ┃ ┃ ┣ 📂cg
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜course_grainer.py
 ┃ ┃ ┃ ┗ 📜multimol_course_grainer.py
 ┃ ┃ ┣ 📂separated
 ┃ ┃ ┃ ┣ 📂parametiser
 ┃ ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┃ ┣ 📜polymer.py
 ┃ ┃ ┃ ┃ ┣ 📜polymer_list.py
 ┃ ┃ ┃ ┃ ┣ 📜solvent.py
 ┃ ┃ ┃ ┃ ┗ 📜solvent_csv.py
 ┃ ┃ ┃ ┗ 📜__init__.py
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┗ 📜base_workflow.py
 ┃ ┣ 📜__init__.py
 ┃ ┗ 📜command_line_operation.py
 ┣ 📂samples
 ┃ ┣ 📜sample_error.csv
 ┃ ┣ 📜sample_output_log.csv
 ┃ ┣ 📜sample_progress.csv
 ┣ 📜README.md
 ┣ 📜main.py
 ┣ 📜parametiser_workflow.py
 ┣ 📜requirements.txt
 ┗ 📜simulation_manager.py
 ```
