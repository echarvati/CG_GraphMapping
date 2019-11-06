#import numpy as np
#import numpy.linalg
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt
from Graphs.MyGraph import MyGraph
import numpy as np

class MyAtom:
    'Collection of the graph nodes'
    def __init__(self, atom, molecule ):
        vertex = molecule.GetAtomWithIdx(atom)
        connectivity = vertex.GetDegree()
        pos = molecule.GetConformer().GetAtomPosition(atom)
        coords = np.array([pos.x, pos.y, pos.z])
        atomtype = molecule.GetAtomWithIdx(atom).GetSymbol()
        # print(type)
        mass = vertex.GetMass()
        chirality = vertex.GetChiralTag()
        charge = vertex.GetFormalCharge()

        self.vertex = vertex
        self.connectivity = connectivity
        self.coords = coords
        self.type = atomtype
        self.weight = mass
        self.chiral = chirality
        self.charge = charge

class Bond:
    'Bonds between atoms or beads'
    def __init__(self, bond_id, bond):
        
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        bd_type = bond.GetBondType()
        if bd_type == Chem.rdchem.BondType.SINGLE:
            bond_type = 1
        elif bd_type == Chem.rdchem.BondType.DOUBLE:
            bond_type = 2
        elif bd_type == Chem.rdchem.BondType.TRIPLE:
            bond_type = 3
        elif bd_type == Chem.rdchem.BondType.AROMATIC:
            bond_type = 1.5
        elif bd_type == None:
            bond_type = None

        self.begin = atom1
        self.end = atom2
        self.id = bond_id
        self.BondOrder = bond_type

class MolecularGraph():
    ''' a wrapper for Graph, atoms, bonds, beads and related function '''
    def __init__(self, SMILES):
        ''' originally mol2graph 
            initialize all thing in CGGraph
            :SMILES: smiles representation of a molecule
        '''
        self.mol_graph = nx.Graph()
        self.others = []
        
        # initialize molecule
        molecule = Chem.MolFromSmiles(SMILES)
        AllChem.EmbedMolecule(molecule, useRandomCoords=True)
        self.molecule = molecule
        
        # initialize atoms, other, bonds and graph
        self.bonds_list = self._build_bonds_graph()
        self.mol_nodes = self._build_atoms()
        self.others = self._build_others()
        self.beads = []

    def _build_atoms(self):
        nodes = []
        for atom in range(len(self.molecule.GetAtoms())):
            nodes.append( MyAtom(atom, self.molecule) )
        return nodes

    def _build_bonds_graph(self):
        bonds_list = []
        for bond_id ,bond in enumerate(self.molecule.GetBonds()):
            this_bond = Bond(bond_id, bond)
            bonds_list.append(this_bond)
            self.mol_graph.add_edge(this_bond.begin, this_bond.end, weight = this_bond.BondOrder)
        return bonds_list

    def _build_others(self):
        others = []
        for item in range(len(self.mol_nodes)):
            # print(mol_nodes[item].type)
            if not self.mol_nodes[item].type == 'C':
                others.append(item)
        return others

    def _single_iteration(self, mygraph):
        ''' this function takes a graph in, analyze, seperate, and contract it
        :graph: big/full MyGraph
        :return: small/contracted networkx graph '''
        # analyze the graph
        vertices, major, minor = mygraph.make_fragment()
        contracted_graph, contracted_nodes = mygraph.contract(minor, major)
        return contracted_graph

    def _make(self):
        mygraph = MyGraph(self.mol_graph, self.others)
        return mygraph

    def iterate(self, n):
        some_graph=self._make()
        for _ in range(n):
            some_graph= MyGraph(self._single_iteration(some_graph))



    