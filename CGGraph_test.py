from Graphs.MolecularGraph import MolecularGraph

# SMILES = 'c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)'
#SMILES = 'c1cc(O)ccc1Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)O'
# SMILES = 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
# SMILES = 'C[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)C[C@H](C)C[C@H](C)'
# SMILES = 'C[C@H](C#N)C[C@@H](C#N)C[C@H](C#N)C[C@@H](C#N)'
SMILES = 'COCCOCCOCCOCCOCCOC'
# SMILES = 'C=CC=CC=CC=CC=CC=C'
# SMILES = 'C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)'
cggraph = MolecularGraph(SMILES)
n = 5
cggraph.iterate(n)