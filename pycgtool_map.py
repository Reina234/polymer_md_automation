from config.data_models.solvent import Solvent


solvent = Solvent("hexane", 86.17, 661, "input_data/solvent_pdbs/hexane.pdb", 1.24e-4)

# from modules.workflows.atomistic.solvent_equilibriator import (
#    SolventEquilibriationWorkflow,
# )

# SolventEquilibriationWorkflow(
#    solvent=solvent, box_size_nm=[7, 7, 7], temperatures=[298], output_dir="1_22_test"
# ).run()

from modules.workflows.atomistic.polymer_equilibriator import (
    PolymerEquilibriationWorkflow,
)


equilibriator = PolymerEquilibriationWorkflow(
    monomer_smiles=["C=Cc1ccccc1", "C=C"],
    num_units=10,
    temperatures=[298],
    solvent=solvent,
    box_size_nm=[7, 7, 7],
    output_dir="ZZZ",
)
outputs = equilibriator.run()

from modules.workflows.cg.course_grainer import CourseGrainer

polymer = equilibriator.polymer

from modules.cg_mappers.pycgtool_map_generator import PyCGToolMapGenerator


mapper = PyCGToolMapGenerator(polymer=polymer)
mapper.add_solvent_to_map(
    itp_file="preprocessed/equilibriated_solvent_boxes/hexane/solvent.itp"
)
mapper.create_map("pycg_map", start_index=1)
