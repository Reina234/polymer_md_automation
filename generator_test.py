from modules.workflows.atomistic.polymer_parametizer import PolymerGeneratorWorkflow

generator = PolymerGeneratorWorkflow(["C=Cc1ccccc1", "C=C"], 10)

generator.run("ZZZ")
