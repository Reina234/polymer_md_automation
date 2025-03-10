# Automated Polymer Simulations

This README is a work in progress.

This repo contains code for automating GROMACS molecular dynamics simulation, currently capable of handling addition homopolymers and simple alternating copolymers.

Overall workflow handling is done via the simulation manager, it uses csv data for solvents, and a list for the monomer smiles (generates all combinations of 1, 2, 3 monomer type composition). This is called in main.py. Caching, error logging, output formatting is all handled by `SimulationManager`. Cache, logs, and temp folders will be automatically generated.

parametiser_workflow.py holds defunct methods that does the whole parameterisation workflow, this was part of inital separation attempts of the parameterisation workflow vs md/gromacs workflow for HPC uploading purposes.

All path inports assumes the repo folder is top level.

## General structure
```
📦config
 ┣ 📂__pycache__
 ┃ ┣ 📜__init__.cpython-312.pyc
 ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┣ 📜acpype_config.cpython-312.pyc
 ┃ ┣ 📜acpype_config.cpython-38.pyc
 ┃ ┣ 📜constants.cpython-312.pyc
 ┃ ┣ 📜constants.cpython-38.pyc
 ┃ ┣ 📜mdp_workflow_config.cpython-38.pyc
 ┃ ┣ 📜paths.cpython-312.pyc
 ┃ ┗ 📜paths.cpython-38.pyc
 ┣ 📂data_models
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜output_types.cpython-38.pyc
 ┃ ┃ ┗ 📜solvent.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜output_types.py
 ┃ ┗ 📜solvent.py
 ┣ 📜__init__.py
 ┣ 📜acpype_config.py
 ┣ 📜constants.py
 ┣ 📜mdp_workflow_config.py
 ┗ 📜paths.py

 📦input_data
 ┣ 📂__pycache__
 ┃ ┣ 📜__init__.cpython-312.pyc
 ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┣ 📜monomer_smiles.cpython-312.pyc
 ┃ ┗ 📜monomer_smiles.cpython-38.pyc
 ┣ 📜__init__.py
 ┣ 📜monomer_smiles.py
 ┣ 📜solvent_data.csv
 ┗ 📜solvent_data.xlsx

 📦modules
 ┣ 📂__pycache__
 ┃ ┣ 📜__init__.cpython-312.pyc
 ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┗ 📜command_line_operation.cpython-38.pyc
 ┣ 📂acpype
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜acpype_parametizer.cpython-38.pyc
 ┃ ┃ ┗ 📜acpype_utils.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜acpype_parametizer.py
 ┃ ┗ 📜acpype_utils.py
 ┣ 📂cache_store
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜base_cache.cpython-38.pyc
 ┃ ┃ ┣ 📜equilibriated_atomistic_polymer_cache.cpython-38.pyc
 ┃ ┃ ┣ 📜file_cache.cpython-38.pyc
 ┃ ┃ ┣ 📜mdp_cache.cpython-38.pyc
 ┃ ┃ ┣ 📜pickle_cache.cpython-38.pyc
 ┃ ┃ ┗ 📜solvent_cache.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜base_cache.py
 ┃ ┣ 📜equilibriated_atomistic_polymer_cache.py
 ┃ ┣ 📜file_cache.py
 ┃ ┣ 📜mdp_cache.py
 ┃ ┣ 📜pickle_cache.py
 ┃ ┗ 📜solvent_cache.py
 ┣ 📂cg_mappers
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜base_map_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜martini_index_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜martini_map_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜multimol_map_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜open_mscg_map_generator.cpython-38.pyc
 ┃ ┃ ┗ 📜pycgtool_map_generator.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜base_map_generator.py
 ┃ ┣ 📜martini_index_generator.py
 ┃ ┣ 📜martini_map_generator.py
 ┃ ┣ 📜multimol_map_generator.py
 ┃ ┣ 📜open_mscg_map_generator.py
 ┃ ┣ 📜pycgtool_map_generator.py
 ┃ ┗ 📜votca_map_generator.py
 ┣ 📂file_conversion
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┗ 📜converter_factory.cpython-38.pyc
 ┃ ┣ 📂converters
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜base_converter.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜editconf_gro_to_pdb.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜editconf_pdb_to_gro.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜obabel_pdb_to_mol2_converter.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜base_converter.py
 ┃ ┃ ┣ 📜editconf_gro_to_pdb.py
 ┃ ┃ ┣ 📜editconf_pdb_to_gro.py
 ┃ ┃ ┗ 📜obabel_pdb_to_mol2_converter.py
 ┃ ┣ 📜__init__.py
 ┃ ┗ 📜converter_factory.py
 ┣ 📂gromacs
 ┃ ┣ 📂equilibriation
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜base_workflow_step.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜full_equilibriation_workflow.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜mdp_cache.cpython-38.pyc
 ┃ ┃ ┣ 📜base_workflow_step.py
 ┃ ┃ ┗ 📜full_equilibriation_workflow.py
 ┃ ┣ 📂parsers
 ┃ ┃ ┣ 📂data_models
 ┃ ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┃ ┗ 📜section.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┗ 📜section.py
 ┃ ┃ ┣ 📂handlers
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜data_handler.py
 ┃ ┃ ┃ ┗ 📜gro_handler.py
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜gromacs_parser.py
 ┃ ┃ ┗ 📜itp_parser.py
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜analyser.py
 ┃ ┗ 📜index_manager.py
 ┣ 📂lammps
 ┃ ┣ 📂parsers
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┗ 📜open_mscg_data_parser.py
 ┃ ┗ 📜__init__.py
 ┣ 📂moltemplate
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┗ 📜moltemplate_utils.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜base_molecule.py
 ┃ ┣ 📜moltemplate_system.py
 ┃ ┣ 📜moltemplate_utils.py
 ┃ ┣ 📜multimol_solvent.py
 ┃ ┣ 📜polymer.py
 ┃ ┗ 📜solvent.py
 ┣ 📂open_mscg
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜force_matcher.cpython-38.pyc
 ┃ ┃ ┣ 📜multimol_topol_exporter.cpython-38.pyc
 ┃ ┃ ┣ 📜multimol_traj_mapper.cpython-38.pyc
 ┃ ┃ ┣ 📜topol_exporter.cpython-38.pyc
 ┃ ┃ ┣ 📜topol_generator.cpython-38.pyc
 ┃ ┃ ┗ 📜trajectory_mapper.cpython-38.pyc
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜force_matcher.py
 ┃ ┣ 📜multimol_topol_exporter.py
 ┃ ┣ 📜multimol_traj_mapper.py
 ┃ ┣ 📜topol_exporter.py
 ┃ ┣ 📜topol_generator.py
 ┃ ┗ 📜trajectory_mapper.py
 ┣ 📂packmol
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜base_packmol_operation.cpython-38.pyc
 ┃ ┃ ┗ 📜solvent_box.cpython-38.pyc
 ┃ ┣ 📂templates
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┗ 📜solvent_box_template.inp
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜base_packmol_operation.py
 ┃ ┗ 📜solvent_box.py
 ┣ 📂rdkit
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📜base_molecule_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜polymer_itp_scaler.cpython-38.pyc
 ┃ ┃ ┗ 📜solvent_generator.cpython-38.pyc
 ┃ ┣ 📂polymer_builders
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜alternating_copolymer.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜base_polymer_generator.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜alternating_copolymer.py
 ┃ ┃ ┣ 📜base_polymer_generator.py
 ┃ ┃ ┗ 📜homopolymer_generator.py
 ┃ ┣ 📜__init__.py
 ┃ ┣ 📜base_molecule_generator.py
 ┃ ┣ 📜polymer_itp_scaler.py
 ┃ ┗ 📜solvent_generator.py
 ┣ 📂utils
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┗ 📜__init__.cpython-38.pyc
 ┃ ┣ 📂atomistic
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜file_utils.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜mdp_utils.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜file_utils.py
 ┃ ┃ ┗ 📜mdp_utils.py
 ┃ ┣ 📂shared
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜calculation_utils.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜dataframe_utils.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜file_utils.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜calculation_utils.py
 ┃ ┃ ┣ 📜dataframe_utils.py
 ┃ ┃ ┗ 📜file_utils.py
 ┃ ┗ 📜__init__.py
 ┣ 📂workflows
 ┃ ┣ 📂__pycache__
 ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┗ 📜base_workflow.cpython-38.pyc
 ┃ ┣ 📂atomistic
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜joined_workflow.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜polymer_equilibriator.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜polymer_parametizer.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜solvent_equilibriator.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜joined_workflow.py
 ┃ ┃ ┣ 📜polymer_equilibriator.py
 ┃ ┃ ┣ 📜polymer_parametizer.py
 ┃ ┃ ┗ 📜solvent_equilibriator.py
 ┃ ┣ 📂cg
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜course_grainer.cpython-38.pyc
 ┃ ┃ ┃ ┗ 📜multimol_course_grainer.cpython-38.pyc
 ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┣ 📜course_grainer.py
 ┃ ┃ ┗ 📜multimol_course_grainer.py
 ┃ ┣ 📂separated
 ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┗ 📜__init__.cpython-38.pyc
 ┃ ┃ ┣ 📂parametiser
 ┃ ┃ ┃ ┣ 📂__pycache__
 ┃ ┃ ┃ ┃ ┣ 📜__init__.cpython-38.pyc
 ┃ ┃ ┃ ┃ ┣ 📜polymer.cpython-38.pyc
 ┃ ┃ ┃ ┃ ┣ 📜polymer_list.cpython-38.pyc
 ┃ ┃ ┃ ┃ ┣ 📜solvent.cpython-38.pyc
 ┃ ┃ ┃ ┃ ┗ 📜solvent_csv.cpython-38.pyc
 ┃ ┃ ┃ ┣ 📜__init__.py
 ┃ ┃ ┃ ┣ 📜polymer.py
 ┃ ┃ ┃ ┣ 📜polymer_list.py
 ┃ ┃ ┃ ┣ 📜solvent.py
 ┃ ┃ ┃ ┗ 📜solvent_csv.py
 ┃ ┃ ┗ 📜__init__.py
 ┃ ┣ 📜__init__.py
 ┃ ┗ 📜base_workflow.py
 ┣ 📜__init__.py
 ┗ 📜command_line_operation.py
 ```