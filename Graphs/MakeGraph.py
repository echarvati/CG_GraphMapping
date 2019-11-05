import numpy as np
import numpy.linalg
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
import networkx as nx


#Define classes
class MyAtom:
    'Collection of the graph nodes'
    def __init__(self, vertex, connect, coordinates, atom_type, mass, chirality, charge ):
        self.vertex = vertex
        self.connectivity = connect
        self.coords = coordinates
        self.type = atom_type
        self.weight = mass
        self.chiral = chirality
        self.charge = charge

class Bond:
    'Bonds between atoms or beads'
    def __init__(self, atom1, atom2, bond_id, order):
        self.begin = atom1
        self.end = atom2
        self.id = bond_id
        self.BondOrder = order

def nodes2vertices(similarity):
    vertices = []
    for sim in range(len(similarity)):
        k = sorted(similarity[sim])
        if k not in vertices:
            vertices.append(k)
        else:
            continue
    sets = []
    common_atoms = []
    for m in range(len(vertices)):
        sets.append(set(vertices[m]))
        for e, sett in enumerate(sets):
            common = set(vertices[m]).intersection(sett)
            if set(vertices[m]) == sett or len(common) == 0:
                continue
            else:
                common_atoms.append(common)
            for c in common:
                if c in sett and c in set(vertices[m]):
                    new = sett.union(set(vertices[m]))
                    sett.update(new)
                    if new.issuperset(set(vertices[m])):
                        sets.remove(set(vertices[m]))
                    else:
                        continue

                else:
                    continue
    vertices.clear()
    for s in sets:
        vertices.append(list(s))

    return vertices

def build_atoms(molecule):
    nodes = []
    for atom in range(len(molecule.GetAtoms())):
        vertex = molecule.GetAtomWithIdx(atom)
        connectivity = vertex.GetDegree()
        pos = molecule.GetConformer().GetAtomPosition(atom)
        coords = np.array([pos.x, pos.y, pos.z])
        type = molecule.GetAtomWithIdx(atom).GetSymbol()
        # print(type)
        mass = vertex.GetMass()
        chirality = vertex.GetChiralTag()
        charge = vertex.GetFormalCharge()
        atom_molecule = MyAtom(vertex, connectivity, coords, type, mass, chirality, charge)
        nodes.append(atom_molecule)
    return nodes

def build_bonds (bond):
    order = []
    atom1 = bond.GetBeginAtomIdx()
    atom2 = bond.GetEndAtomIdx()
    type = bond.GetBondType()
    if type == Chem.rdchem.BondType.SINGLE:
        bond_type = 1
    elif type == Chem.rdchem.BondType.DOUBLE:
        bond_type = 2
    elif type == Chem.rdchem.BondType.TRIPLE:
        bond_type = 3
    elif type == Chem.rdchem.BondType.AROMATIC:
        bond_type = 1.5
    elif type == None:
        pass

    return atom1, atom2, bond_type


def find_rings(graph, rings):
    temp = []
    backbone = []
    roots = []

    for ring in rings:
        for ring_nd in ring:
            temp.append(ring_nd)
    for node in graph.nodes():
        if node in temp:
            continue
        else:
            backbone.append(node)
    for node_bck in backbone:
        for node_r in temp:
            if graph.has_edge(node_bck, node_r):
                roots.append(node_bck)
            else:
                continue

    return rings, roots, backbone

def find_branches(graph):
    graph_branches = []
    test_brnh = []
    temp_branches = []
    roots = []
    start = list(graph.nodes())[0]
    stop = len(graph) - 1

    backs = list(nx.shortest_simple_paths(graph, source=start, target=stop))
    for back in backs:
        backbone = back

    for node_gr in graph.nodes():
        if node_gr in backbone:
            continue
        else:
            temp_branches.append(node_gr)

    for back_nd in backbone:
        for brcnh_nd in temp_branches:
            if graph.has_edge(back_nd, brcnh_nd):
                roots.append(back_nd)
            else:
                continue
        for num, temp_brcnd in enumerate(temp_branches):
            nd1 = temp_brcnd
            nd2 = temp_branches[(num+1)%len(temp_branches)]
            if graph.has_edge(nd1, nd2):
                    test_brnh.append([nd1, nd2])
                    graph_branches = nodes2vertices(test_brnh)
            else:
                for root in roots:
                    if graph.has_edge(nd1, root) and graph.degree(nd1) == 1 and [nd1] not in graph_branches:
                        graph_branches.append([nd1])

    return graph_branches, roots, backbone

def find_regions(graph):

    matches_rings = nx.cycle_basis(graph)
    if len(matches_rings)!=0:
       rings, roots, backbone = find_rings(graph, matches_rings)
       branches = []
    else:
        branches, roots, backbone = find_branches(graph)
        rings = []

    return rings, branches, roots, backbone

def find_rigid(graph, region):
    conj_func = []

    for (u, v, d) in graph.edges(data=True):
        if u in region and v in region:
            if graph.has_edge(u, v):

                if graph[u][v]['weight'] == 1:
                    conj_func.append(-1)
                elif graph[u][v]['weight'] > 1:
                    conj_func.append(1)

    return conj_func

def conj_region(conj):

    conj_func = sum(conj)
    if conj_func>=0:
        rigid=True
    else:
        rigid=False
    return rigid


def check_rigid(graph, rings, branches, backbone):
    rigid_ring = 0
    rigid_branch = 0
    rigid_backbone = 0

    if len(rings)!=0:
        for ring in rings:
            ring_conj_func = find_rigid(graph, ring)
        rigid_ring = conj_region(ring_conj_func)
    elif len(branches)!=0:
        for branch in branches:
            branch_conj_func = find_rigid(graph, branch)
        rigid_branch = conj_region(branch_conj_func)
    else:
        backbone_conj_func = find_rigid(graph, backbone)
        rigid_backbone = conj_region(backbone_conj_func)

    return rigid_ring, rigid_branch, rigid_backbone

def mol2graph(SMILES):
    mol_graph = nx.Graph()
    mol_nodes = []
    bonds_list = []
    others = []

    molecule = Chem.MolFromSmiles(SMILES)
    AllChem.EmbedMolecule(molecule, useRandomCoords=True)
    mol_nodes = build_atoms(molecule)
    for item in range(len(mol_nodes)):
        # print(mol_nodes[item].type)
        if mol_nodes[item].type == 'C':
            continue
        else:
            others.append(item)

    for num,bond in enumerate(molecule.GetBonds()):
        bond_id = num
        atom1, atom2, order = build_bonds(bond)
        molec_bond = Bond(atom1, atom2, bond_id, order)
        bonds_list.append(molec_bond)
        mol_graph.add_edge(atom1, atom2, weight = order)

    rings, branches, roots, backbone = find_regions(mol_graph)

    rigid_ring, rigid_branch, rigid_backbone = check_rigid(mol_graph, rings, branches, backbone)

    return mol_graph, mol_nodes, rings, branches, roots, backbone, rigid_ring, rigid_branch, rigid_backbone, others
