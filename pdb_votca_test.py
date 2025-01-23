from rdkit_new.votca_pdb_generator import VOTCASingleBeadMapping
from data_models.solvent import Solvent

solvent = Solvent("hexane", 86.17, 661, "input_data/solvent_pdbs/hexane.pdb", 1.24e-4)

generator = VOTCASingleBeadMapping(
    molecule_name=solvent.pdb_molecule_name,
    itp_file="1_22_test/hexane/solvent.itp",
)
generator.generate_xml("solvent.xml")
