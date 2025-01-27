from modules.rdkit.polymer_builders.alternating_copolymer import (
    AlternatingPolymerGenerator,
)
from modules.rdkit.polymer_itp_scaler import PolymerITPScaler

generator = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator.generate_polymer(10, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

generator2 = AlternatingPolymerGenerator(["C=Cc1ccccc1"])
pdb = generator2.generate_polymer(3, "rdkit_test3", overwrite=True, save=True)
# parameterize(pdb, "alternating_test")

itp = PolymerITPScaler(
    "1_22_test/c=cc1ccccc1_3/POLY_GMX.itp",
    generator2.cg_map,
    7,
)

itp.write_itp("test_3.itp")
