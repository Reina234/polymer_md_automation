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

from modules.workflows.cg.multimol_course_grainer import MultimolCourseGrainer

polymer = equilibriator.polymer
print(polymer.cg_map)
from modules.workflows.cg.course_grainer import CourseGrainer

# cg = CourseGrainer(polymer=equilibriator.polymer, outputs=outputs[0], output_dir="ZZZ2")
# cg.run()

# cg2 = MultimolCourseGrainer(
#    polymer=equilibriator.polymer, outputs=outputs[0], output_dir="ZZZ"
# )
# cg2.run()


from modules.moltemplate.moltemplate_system import MoltemplateSystem

system = MoltemplateSystem(
    n_units=100,
    polymer=polymer,
    box_nm=[7, 7, 7],
    solvent=solvent,
    openmscg_topol_path="cg_poly.data",
)

system.write_system_lt()
