from rdkit import Chem

ptable = Chem.GetPeriodicTable()
for i in range(1, ptable.GetElements() + 1):
    element = Chem.GetPeriodicTable().GetElementSymbol(i)
    print
all_elements = [ptable.GetElementSymbol(i) for i in range(1, ptable.GetElements() + 1)]

print("Recognized elements in RDKit:", all_elements)
