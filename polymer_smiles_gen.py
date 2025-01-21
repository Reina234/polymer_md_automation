from modules.atomistic.rdkit.homopolymer_generator import HomopolymerGenerator

generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 5, "rdkit_test2", overwrite=False, save=False)
print(generator.end_caps)
