from A_modules.atomistic.rdkit.homopolymer_generator import HomopolymerGenerator

generator = HomopolymerGenerator()
generator.generate_polymer("C=Cc1ccccc1", 5, "rdkit_test")
