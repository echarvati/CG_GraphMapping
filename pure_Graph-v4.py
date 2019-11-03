#Graph
import numpy as np
import numpy.linalg
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt
import time
from Graphs.DrawGraph import draw_graph as draw_graph
from Graphs.MakeGraph import mol2graph as mol2graph
from Graphs.MakeFragments import make_fragment as mapping
from Graphs.MakeFragments import contract as contract

start_time = time.clock()



# SMILES = 'C=CC=CC=CC=CC=CC=CC'
SMILES = 'CCCCCCCCCCCCCCC'
# SMILES = 'COCOCOCOCOCOC'
# SMILES = 'c1cc(O)ccc1Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)O'
# SMILES = 'C[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)C[C@H](C)C'
# SMILES = 'C[C@H](C#N)C[C@@H](C#N)C[C@H]C'
# SMILES = 'c1ccccc1[C@H](C)C[C@@H](c1ccccc1)C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C'
# SMILES = 'c1ccccc1[C@H](C)C[C@@H](c1ccccc1)C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C'
#SMILES = 'C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C'
# SMILES = 'C(F)C(F)C(F)C(F)C(F)C(F)C(F)'
#  SMILES = 'CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)C'
# SMILES = 'C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)'
#SMILES = 'c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)c1ccc(cc1)'
#SMILES = 'c1cc(O)ccc1Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)Oc1ccc(cc1)O'
# SMILES = 'c1ccccc1[C@H](C)C[C@@H](c1ccccc1)C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C'

AAgraph, AAnodes, rings, branches, roots, backbone, rigid_ring, rigid_branch, rigid_backbone, other_atoms= mol2graph(SMILES)


vertices, major, minor = mapping(AAgraph, AAnodes, rings, branches, roots, backbone, rigid_ring, rigid_branch, rigid_backbone, other_atoms)
# print(vertices)
# print(major)
# print(minor)
draw_graph(AAgraph, AAnodes, other_atoms)
CGgraph, CGnodes = contract(minor, major, AAgraph)
print(CGnodes)
draw_graph(CGgraph)


#