import numpy as np
import numpy.linalg
from rdkit import Chem
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt
import time

start_time = time.clock()

#Graph
import numpy as np
import numpy.linalg
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
import networkx as nx
import matplotlib.pyplot as plt
import time

start_time = time.clock()

#Graph related operations

#mol2graph: Read the SMILES string input as molecule (RDKit), extract all relevant information (mass, coordinates,
#chirality, charge etc) and save it on MyAtom class. The molecule is, then, converted to a molecular graph graph
#where the atoms are the graph nodes, the bonds the graph edges and the bond order the edge weight. Information such as
#chirality and charges is does not yet appear on the  molecular graph.

def mol2graph(SMILES):
    mol_graph = nx.Graph()
    mol_nodes = []

    nodes_list = []
    edges_list = []
    bonds_list = []

    molecule = Chem.MolFromSmiles(SMILES)
    # mol = Chem.AddHs(molecule)
    AllChem.EmbedMolecule(molecule, useRandomCoords=True)

    for atom in range(len(molecule.GetAtoms())):
        pos = molecule.GetConformer().GetAtomPosition(atom)
        coords = np.array([pos.x, pos.y, pos.z])
        vertex = molecule.GetAtomWithIdx(atom)
        neigh = vertex.GetDegree()
        mass = vertex.GetMass()
        type = vertex.GetSymbol()
        chirality = vertex.GetChiralTag()
        charge = vertex.GetFormalCharge()
        atom_molecule = MyAtom(vertex, coords, type, mass, chirality, charge, neigh)
        mol_nodes.append(atom_molecule)


    order = []

    for bond in molecule.GetBonds():
        atom1 = bond.GetBeginAtomIdx()
        atom2 = bond.GetEndAtomIdx()
        type = bond.GetBondType()
        if type == Chem.rdchem.BondType.SINGLE:
            bond_type = 1
            order.append(1)
        elif type == Chem.rdchem.BondType.DOUBLE:
            bond_type = 2
            order.append(2)
        elif type == Chem.rdchem.BondType.TRIPLE:
            bond_type = 3
            order.append(3)
        elif type == Chem.rdchem.BondType.AROMATIC:
            bond_type = 1.5
            order.append(1.5)
        elif type == None:
            continue
        bond_molec = Bond(atom1, atom2, bond_type)
        bonds_list.append(bond_molec)
        mol_graph.add_edge(atom1, atom2)

    for e1, e2 in mol_graph.edges():
        o = order[e1]
        mol_graph[e1][e2]['weight'] = o

    w = " ".join(str(e) for e in order)

    carbons = []
    others = []

    for a in bonds_list:
        type_1 = mol_nodes[a.begin].type
        if type_1 == 'C':
            carbons.append(a.begin)
        else:
            others.append(a.begin)

        type_2 = mol_nodes[a.end].type
        if type_2 == 'C' and a.end != carbons:
            carbons.append(a.end)
        else:
            if type_2 != 'C' and a.end not in others:
                others.append(a.end)

    carbons = list(dict.fromkeys(carbons))
    others = list(dict.fromkeys(others))


    return(mol_graph, mol_nodes, carbons, others)


def draw_graph(graph, nodes=None, atom_type1=None, atom_type2=None):
    color_map = []

    esingle = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 1]
    etriple = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 3]
    edouble = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 2]
    earom = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 1.5]
    # var = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 0.75]

    pos = nx.spring_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=500)
    # nx.draw_networkx_nodes(graph, pos, node_size=250, nodelist=atom_type1, node_color='g')
    if atom_type2 != None and nodes!=None:
        for atom in atom_type2:
            if nodes[atom].type == 'O':
                nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=atom_type2, node_color='r')
            elif nodes[atom].type == 'Cl':
                nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=atom_type2, node_color='y')
            elif nodes[atom].type == 'F':
                nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=atom_type2, node_color='g')

    nx.draw_networkx_edges(graph, pos, edgelist=etriple, width=12, alpha=0.5, edge_color='b', node_color=color_map)
    nx.draw_networkx_edges(graph, pos, edgelist=edouble, width=12, alpha=0.5, edge_color='g', node_color=color_map)
    nx.draw_networkx_edges(graph, pos, edgelist=esingle, width=10, node_color=color_map)
    nx.draw_networkx_edges(graph, pos, edgelist=earom, width=12, alpha=0.5, edge_color='r', node_color=color_map)
    # nx.draw_networkx_edges(graph, pos, edgelist=var, width=3, style='dashed')
    nx.draw_networkx_labels(graph, pos, font_size=20, font_family='sans-serif', node_color=color_map)
    plt.axis('off')
    plt.show()

def check_features(graph, atom_type1=None):
    cycle = []
    matches_rings = nx.cycle_basis(graph)
    neighbors = []
    temp = []
    back = []
    backbone = []
    branches = []
    roots = []
    leaf=[]

    if len(matches_rings) != 0:
        for match in matches_rings:
            for elm in match:
                cycle.append(elm)
        for n in graph.nodes():
            neighbors.append(list(graph.neighbors(n)))
            for neighbor in neighbors:
                for n in neighbor:
                    temp.append(n)
                    if n in cycle:
                        temp.remove(n)
                    else:
                        continue
        temp = list(dict.fromkeys(temp))
        # print(temp)
        back.append(temp)
        backbone = back
        for bck in back:
            for b in bck:
                for ccl in cycle:
                    if graph.has_edge(b,ccl):
                        roots.append(b)
    else:
        start = list(graph.nodes())[0]
        stop = len(graph) - 1
        backbone = list(nx.shortest_simple_paths(graph, source=start, target=stop))
        for node in list(graph.nodes()):
            for babone in backbone:
                if node in babone:
                    continue
                else:
                    for bone in backbone:
                        for nd in graph.nodes():
                            if nd not in bone:
                                for bon in bone:
                                    if graph.has_edge(bon, nd) == True:
                                        roots.append(bon)
                                        if graph.degree[nd] == 1:
                                            if [nd] not in branches:
                                                branches.append([nd])

                                        else:
                                            T = nx.dfs_tree(graph)
                                            smth = T.edges()
                                            for s in smth:
                                                a, b = s
                                                if a not in bone and b not in bone:
                                                        if graph.has_edge(a, b) == True:
                                                            leaf.append([a, b])
                                                            branches = nodes2vertices(leaf)
    for back in backbone:
        for a in back:

            if atom_type1!=None and back[-1] not in atom_type1:
                if graph.has_edge(a, back[-1]):
                    roots.append(a)
                else:
                    continue
                branches.append([back[-1]])
                back.remove(back[-1])
            else:
                continue

    roots = list(dict.fromkeys(roots))
    return matches_rings, branches, roots, backbone

def check_branch(branches):
    for branch in branches:
        if len(branch)>1:
            big_branch = True
        else:
            big_branch = False
    return big_branch


def check_rigid(graph, feature):
    global rigid
    conj_func = []
    for (u, v, d) in graph.edges(data=True):
        for feat in feature:
            if u in feat and v in feat:
                if graph.has_edge(u, v) == True:
                    if graph[u][v]['weight'] != 1:
                        conj_func.append(1)
                    else:
                        conj_func.append(-1)

    conj = sum(conj_func)
    if conj>=0:
        rigid=True
    else:
        rigid=False

    return rigid

def make_SDKvertex(graph, atom_type1, atom_type2, backbone, roots, branches):
    if len(atom_type2) == 0:
        if len(branches) == 0:
            vertex = SDK_simple(backbone)
            return vertex
        else:
            backbone_vertex = SDK_branch(backbone, roots, graph)
            branch_vertex = SDK_simple(branches)
            return backbone_vertex, branch_vertex
    else:
        groups = SDK_heavy(graph, atom_type1, atom_type2, backbone)
        return groups

def RigidBond_similarity(graph, in_list):
    sim_nd = []
    diff_nd = []
    major_nd = []
    minor_nd = []

    Lap = nx.normalized_laplacian_matrix(graph, nodelist=None)
    print(Lap)
    for node in in_list:
        for nd in node:
            for oe in node:
                if nd == oe:
                    continue
                else:
                    sim_rank = abs(Lap.A[nd][oe])
                    # print(sim_rank, nd, oe)
                    if sim_rank == 0 or sim_rank == 1:
                        continue
                    elif sim_rank > 0.5:
                        sim_nd.append([nd, oe])
                    elif sim_rank < 0.5:
                        diff_nd.append([nd, oe])
    vertices = nodes2vertices(sim_nd)

    pagerank = nx.pagerank_numpy(graph, alpha=0.85, weight='w')
    print(pagerank)
    (print(nx.eigenvector_centrality_numpy(graph)))

    for bead in vertices:
        mj, mn, tch = node_rank(bead, pagerank, graph)
        for m in mj:
            major_nd.append(m)
        for n in mn:
            minor_nd.append(n)

    return vertices, major_nd, minor_nd

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

def SDK_simple(backbone, graph):
    chn = []
    sub_chn = []
    groups = []
    for chain in backbone:
        if len(chain) <= 3:
            groups.append(chain)
        elif len(chain) == 4:
            for c in range(0, len(chain), 2):
                chn = chain[c: c + 2]
                groups.append(chn)
        elif len(chain) > 4:
            if len(chain) % 3 == 0:
                for c in range(0, len(chain), 3):
                    chn = chain[c: c + 3]
                    groups.append(chn)
            else:
                sub_chn.append([chain[0], chain[1]])
                # for sb in sub_chn:
                #     groups.append(sb)
                #     # for s in sb:
                    #     chain.remove(s)
                for c in range(0, len(chain), 3):
                    chn = chain[c: c + 3]
                    if chn not in groups:
                        groups.append(chn)
    major_nd = []
    minor_nd = []

    for grp in groups:
        sub_graph = nx.Graph(graph.subgraph(grp))
        SUB_pagerank = nx.pagerank_numpy(sub_graph, alpha=0.85, weight='w')
        print(grp)
        print(nx.betweenness_centrality(sub_graph, weight='w'))
        print(nx.eigenvector_centrality(sub_graph))
        print(sub_graph.degree(weight='w'))
        print(nx.adjacency_matrix(sub_graph))
        print(nx.normalized_laplacian_matrix(sub_graph))
        print(SUB_pagerank)
        mj, mn, tch = node_rank(grp, SUB_pagerank, sub_graph)
        for m in mj:
            major_nd.append(m)
        for n in mn:
            minor_nd.append(n)

    return groups, major_nd, minor_nd

def SDK_branch(backbone, branches, roots, graph):
    groups_back = []
    major_nd = []
    minor_nd = []

    if len(branches)!=0:
        for back in backbone:
            if roots == back:
                groups_back = SDK_rooted(branches, backbone, back)
            else:
                if check_branch(branches)==False:
                    groups_back = SDK_small(graph, roots)
                else:
                    groups_back = SDK_tree(graph, backbone, roots)
    else:
        groups_back = SDK_tree(graph, backbone, roots)

    for grp in groups_back:
        sub_graph = nx.Graph(graph.subgraph(grp))
        print(grp)
        print(nx.betweenness_centrality(sub_graph, weight='w'))
        print(nx.eigenvector_centrality(sub_graph))
        print(sub_graph.degree(weight = 'w'))
        print(nx.adjacency_matrix(sub_graph))
        print(nx.normalized_laplacian_matrix(sub_graph))

        SUB_pagerank = nx.pagerank_numpy(sub_graph, alpha=0.85, weight='w')
        print(SUB_pagerank)

        mj, mn, tch = node_rank(grp, SUB_pagerank, sub_graph)
        for m in mj:
            major_nd.append(m)
        for n in mn:
            minor_nd.append(n)
    return groups_back, major_nd, minor_nd

def SDK_rooted(branches, graph, back):
    groups = []
    smth = []
    for branch in branches:
        for at1 in branch:
            smth.append(at1)
    print(smth)

    for at in back:
        for s1 in smth:
            for s2 in smth:
                if s1!=s2:
                    if graph.has_edge(at, s1) and graph.has_edge(at, s2) and sorted([at, s1, s2]) not in groups:
                        groups.append(sorted([at, s1, s2]))
                    elif graph.has_edge(at, s1) and sorted([at, s1]) not in groups:
                        groups.append(sorted([at, s1]))


    return groups


def SDK_small(graph, roots):
    groups = []
    for root in roots:
        small = list(nx.neighbors(graph, root))
        small.append(root)
        if sorted(small) not in groups:
            groups.append(sorted(small))

    return groups


def SDK_tree(graph, backbone, roots):
    short_back = []
    groups = []

    for smth in backbone:
        for s in smth:
            if s in roots:
                continue
            else:
                short_back.append(s)
    # print(short_back)

    for bc, ak in enumerate(short_back):
        p1 = bc
        p2 = short_back[(bc + 1) % len(short_back)]
        if short_back[p1] == short_back[-1]:
            continue
        else:
            path = nx.shortest_path(graph, source=short_back[p1], target=p2)
            if len(path) < 3:
                continue
            elif path not in groups:
                groups.append(path)

    # print(groups)
    return groups

def SDK_heavy(graph, atom_type1, atom_type2, backbone):
    tch = []
    sub_grps = []
    groups = []
    for chain in backbone:
        for c in atom_type1:
            for c_t, c_n in enumerate(atom_type2):
                ths = c_t
                nxt = atom_type2[(c_t + 1) % len(atom_type2)]
                if graph.has_edge(atom_type1[1], atom_type2[ths]) == True and [atom_type1[0], atom_type1[1],
                                                                               atom_type2[0]] not in sub_grps:
                    sub_grps.append([atom_type1[0], atom_type1[1], atom_type2[0]])
                elif atom_type2[ths] == atom_type2[-1] and graph.has_edge(atom_type2[-1], c) == True \
                        and [atom_type2[-1], atom_type1[-1]] not in sub_grps:
                    sub_grps.append(sorted(nx.shortest_path(graph, source=atom_type2[-1], target=atom_type1[-1])))
                elif atom_type2[ths] != atom_type2[-1] and [atom_type2[ths], nxt] not in tch:
                    sub_grps.append(nx.shortest_path(graph, source=atom_type2[ths], target=nxt))
                    tch.append([atom_type2[ths], nxt])

    from operator import itemgetter
    from itertools import groupby

    sub_grps.sort(key=itemgetter(1))

    for elt, items in groupby(sub_grps, itemgetter(1)):
        for i in items:
            groups.append(i)

    major_nd = []
    minor_nd = []

    for grp in sub_grps:
        sub_graph = nx.Graph(graph.subgraph(grp))
        SUB_pagerank = nx.pagerank_numpy(sub_graph, alpha=0.85, weight='w')
        print(grp)
        print('sub degree', sub_graph.degree)
        print('fragment eigenvector', nx.eigenvector_centrality(sub_graph))
        print('fragment betweenness', nx.betweenness_centrality(sub_graph))
        print('Adj', nx.adjacency_matrix(sub_graph))
        print('Lap', nx.normalized_laplacian_matrix(sub_graph))

        mj, mn, tch = node_rank(grp, SUB_pagerank, sub_graph)
        for m in mj:
            major_nd.append(m)
        for n in mn:
            minor_nd.append(n)

    return groups, major_nd, minor_nd

def position_heavy(atom_type2, backbone, branches, rings):
    for at2 in atom_type2:
        for back in backbone:
            for branch in branches:
                for ring in rings:
                    if at2 in back:
                        heavy_back = True
                        heavy_branch = False
                        heavy_ring = False
                    elif at2 in branch:
                        heavy_back = False
                        heavy_branch = True
                        heavy_ring = False
                    elif at2 in ring:
                        heavy_back = False
                        heavy_branch = False
                        heavy_ring = True

    return heavy_back, heavy_branch, heavy_ring


def make_SDKvertex(graph, atom_type1, atom_type2, backbone, roots, branches, rings):

    for babone in backbone:
        for at in babone:
            if at in atom_type2:
                # print(at)
                groups, major_nd, minor_nd = SDK_heavy(graph, atom_type1, atom_type2, backbone)
                return groups, major_nd, minor_nd
            else:
                continue

    if len(branches) == 0 and len(rings)==0:
        groups, major_nd, minor_nd  = SDK_simple(backbone, graph)
        return groups, major_nd, minor_nd
    elif len(branches)!=0 and len(rings)==0:
        for branch in branches:
            if len(branch)==1:
                # print('branch')
                backbone_vertex, major_vertex, minor_vertex = SDK_branch(backbone, branches, roots, graph)
                return backbone_vertex, major_vertex, minor_vertex
            elif len(branch)>1:
                # print('branch')
                backbone_vertex, major_vertex, minor_vertex = SDK_branch(backbone, branches, roots, graph)
                branch_vertex, major_branch, minor_branch = SDK_simple(branches, graph)
                return backbone_vertex, major_vertex, minor_vertex, branch_vertex, major_branch, minor_branch
    elif len(branches)==0 and len(rings)!=0:
        backbone_vertex, major_vertex, minor_vertex = SDK_branch(backbone, branches, roots, graph)
        return backbone_vertex, major_vertex, minor_vertex

def node_rank (grp, centrality, graph):
    major = []
    minor = []
    touched = []
    for g in grp:
        for p in grp:
            if g==p:
                continue
            else:
                if g not in touched and p not in touched and graph.has_edge(g,p):
                    if round(centrality[g], 3) == round(centrality[p], 3):
                        touched.append(g)
                        touched.append(p)
                    elif round(centrality[g], 3) > round(centrality[p], 3):
                        major.append(g)
                        minor.append(p)
                    elif round(centrality[g], 3) < round(centrality[p], 3):
                        major.append(p)
                        minor.append(g)
                else:
                    continue
            for m in major:
                for j in major:
                    if m in grp and j in grp and m!=j:
                        major.remove(j)
                        minor.append(j)
                    else:
                        continue
        if g not in major and g not in minor and p not in major and g not in minor:
            if g!=p:
                major.append(g)
                minor.append(p)
        else:
            continue
    major = list((dict.fromkeys(major)))


    return(major, minor, touched)


def contract_minor_nodes(minor, graph):
    CGgraph = graph.copy()
    for mn_node in minor:
        if mn_node in CGgraph.nodes():
            CGgraph.remove_node(mn_node)
        else:
            continue
    return CGgraph

#Make beads and update nodes
def make_beadmass(match, shared, nodes):
    m_atoms = []
    p_atoms = []
    atoms = []
    for m in match:
        if m in shared:
            m_atoms.append(nodes[m].weight * 0.5)
            p_atoms.append(nodes[m].coords * 0.5)
        else:
            m_atoms.append(nodes[m].weight)
            p_atoms.append(nodes[m].coords)
        atoms.append(nodes[m].content)
        mtot = sum(m_atoms)

    return atoms, p_atoms, m_atoms, mtot

def make_beadcoords(p_atoms, m_atoms, mbead):
    x = 0.0
    y = 0.0
    z = 0.0
    for k in range(len(p_atoms)):
        x = x + p_atoms[k][0] * m_atoms[k]
        y = y + p_atoms[k][1] * m_atoms[k]
        z = z + p_atoms[k][2] * m_atoms[k]
    x = float(x) / mbead
    y = float(y) / mbead
    z = float(z) / mbead
    CoM = np.array([x, y, z])

    return CoM

def make_beadid(match, nodes, atom_type1, atom_type2, backbone, graph, major):
    counter1 = 0
    counter2 = 0
    types = []
    bead_position = []
    position_tag = ''
    chiral_tag = ''
    type1= ''
    typed = [ ]
    for m in match:
        print(m, match, nodes[m].type)
        type = nodes[m].type
        types.append(type)
        for bck in backbone:
            if m in bck:
                if m == bck[0] or m==bck[-1]:
                    bead_position.append('E')
                elif graph.degree(m)>1:
                    bead_position.append('M')
            else:
                continue
        if m in major:
            if hasattr(nodes[m], 'chiral'):
                if nodes[m].chiral == rdkit.Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                    chiral_tag = 'A'
                elif nodes[m].chiral == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW:
                    chiral_tag = 'S'
                elif nodes[m].chiral == rdkit.Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW:
                    chiral_tag = 'R'
            else:
                pass

    if(len(set(types))==1):
        type1 = types[0]
        counter1 = counter1 + len(types)
    else:
        pass

    if len(bead_position)!=0:
        if 'E' in bead_position:
            position_tag = 'E'
        else:
            position_tag = 'M'
    else:
        pass

    if len(atom_type2)!=0:
        name = type1, str(counter1), type2, str(counter2), position_tag, chiral_tag
        id = "".join(name)
    else:
        name =  type1, str(counter1), position_tag, chiral_tag
        id = "".join(name)

    return(id)


def make_bead(matches, nodes, atom_type1, atom_type2, backbone, graph, major):
    beads_list = []
    shared_atoms = []

    for g1, g2 in enumerate(matches):
        ths_grp = g1
        nxt_grp = matches[(ths_grp + 1) % len(matches)]
        for g in matches[ths_grp]:
            if g in matches[ths_grp] and g in nxt_grp:
                shared_atoms.append(g)
    for match in matches:
        bead_id = make_beadid(match, nodes, atom_type1, atom_type2, backbone, graph, major)
        print(match, bead_id)
        atoms, p_atoms, m_atoms, mtot= make_beadmass(match, shared_atoms, nodes)
        bead_CoM = make_beadcoords(p_atoms, m_atoms, mtot)

        group = Bead(match, bead_CoM, mtot, bead_id)
        beads_list.append(group)

    return beads_list

def update_nodes_aftergraining(bead_matches, CGgraph, beads_list, nodes):
    CG_nodes = nodes.copy()
    for match in bead_matches:
        for am in match:
            if am in CGgraph.nodes():
                CG_nodes[am] = 'BeadPosition'
            else:
                CG_nodes[am] = 'NullVoid'
    tchd = []
    for b in range(len(beads_list)):
        for n in range(len(CG_nodes)):
            if CG_nodes[n] == 'BeadPosition':
                tchd.append(beads_list[b])
                CG_nodes[n] = beads_list[b]
                if beads_list[b] in tchd:
                    break
                else:
                    continue
    while 'NullVoid' in CG_nodes:
        CG_nodes.remove('NullVoid')
    while 'BeadPosition' in CG_nodes:
        CG_nodes.remove('BeadPosition')
    return CG_nodes

def make_xyz(filename, nodes):
    out = open(filename, 'w')
    part_tot = len(nodes)

    print(part_tot, file=out)
    print('CG contraction', SMILES, file=out)

    for i in range(len(nodes)):
        k = " ".join(repr(e) for e in nodes[i].coords)
        part_type, part_coords = (nodes[i].type, k)
        print(part_type, part_coords, file=out)


#Define classes
class MyAtom:
    'Collection of the graph nodes'
    def __init__(self, vertex, coordinates, atom_type, mass, chirality, charge, connect):
        self.content = vertex
        # self.centrality = centrality
        self.coords = coordinates
        self.type = atom_type
        self.weight = mass
        self.chiral = chirality
        self.charge = charge
        self.connectivity = connect

class Bond:
    'Bonds between atoms or beads'
    def __init__(self, atom1, atom2,  order):
        self.begin = atom1
        self.end = atom2
        self.BondOrder = order

class Bead:
#     'Collection of double bonded beads'
    def __init__(self, bead_nodes, bead_CoM,  bead_mass, bead_id):
        self.content = bead_nodes
        self.coords = bead_CoM
        self.weight = bead_mass
        self.type = bead_id


##########################Main
# SMILES = 'C#CC[C@H](C[C@H](/C=C/C)c1ccccc1)c1ccccc1'
# SMILES = 'c1ccccc1[C@H](C)C[C@@H](c1ccccc1)C[C@H](c1ccccc1)C[C@@H](c1ccccc1)C'
# SMILES = 'C[C@H](C)C[C@@H](C)C[C@H](C)C[C@@H](C)C[C@H](C)C'
# SMILES = 'C[C@H](Cl)C[C@@H](Cl)C[C@H](Cl)C[C@@H](Cl)C[C@H](Cl)C'
# SMILES = 'C[C@H](C#N)C[C@@H](C#N)C[C@H](C#N)C[C@@H](C#N)C[C@H](C#N)C'
# SMILES = 'C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)C(F)(F)'
# SMILES = 'C(F)C(F)C(F)C(F)C(F)C(F)C(F)'
# SMILES = 'C(C)C(C)C(C)C(C)C(C)'
# SMILES = 'C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C'
SMILES = 'C[C@H](C#N)C[C@H](C#N)C[C@H](C#N)C[C@@H](C#N)C[C@H](C#N)C[C@H]'
# SMILES = 'CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)CC(C)(CCC)C'
# SMILES = 'CCCCCCCCCCCC'
# SMILES = 'CCOCOCOCOCOCOCOC'
# SMILES = 'CCOCCOCCOCCOCCOCCOCCOCC'
# SMILES = 'C=CC=CC=CC=CC=CC=CC=CC=CC=CC=C'
# SMILES = 'C=CC=CC=CC=CC=CC=CC=C'
# SMILES = 'C=CC=CC=CC=CC=CC=CC'
# SMILES = 'C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C[C@@H](CCC)C[C@H](CCC)C'

AA_molecule, AA_nodes, atom_type1, atom_type2 = mol2graph(SMILES)

draw_graph(AA_molecule, AA_nodes, atom_type1, atom_type2)

AArings, AAbranches, AAroots, AAbackbone = check_features(AA_molecule, atom_type1)

print('carbons', atom_type1)
print('fluorine', atom_type2)
print('rings', AArings)
print('backbone', AAbackbone)
print('branches', AAbranches)
print('roots', AAroots)

bond_backbone = check_rigid(AA_molecule, AAbackbone)
# print(check_rigid(AA_molecule, AA_backbone))

if len(AArings)!=0:
    bond_ring = check_rigid(AA_molecule, AArings)
    # print(bond_ring)
if len(AAbranches)!=0:
    big_branch = check_branch(AAbranches)
    if big_branch==True:
        bond_branch = check_rigid(AA_molecule, AAbranches)
    else:
        pass

# print('big branches?', big_branch)
if bond_backbone == True:
    similar_back, major_back, minor_back = RigidBond_similarity(AA_molecule, AAbackbone)
    major = major_back
    minor = minor_back
    vertices = similar_back
else:
    if len(AAbranches)!=0:
        if big_branch == False:
            print('SDK')
            babr_vertex, babr_major, babr_minor = make_SDKvertex(AA_molecule, atom_type1, atom_type2, AAbackbone, AAroots, AAbranches, AArings)
            print(babr_vertex)
            print(babr_major)
            print(babr_minor)
            major = babr_major
            minor = babr_minor
            vertices = babr_vertex

        elif bond_branch == False and big_branch == True:
            backbone_vertex, backbone_major, backbone_minor, branch_vertex, branch_major, branch_minor = make_SDKvertex(AA_molecule, atom_type1, atom_type2, AAbackbone, AAroots, AAbranches, AArings)
            print('backbone', backbone_vertex)
            print('major back', backbone_major)
            print('minor back', backbone_minor)
            print(branch_vertex)
            print(branch_major)
            print(branch_minor)
            major = backbone_major + branch_major
            minor = backbone_minor + branch_minor
            vertices = backbone_vertex + branch_vertex

    elif len(AArings)!=0:
        if bond_backbone==False and bond_ring==True:
            similar_ring, major_ring, minor_ring = RigidBond_similarity(AA_molecule, AArings)
            backbone_vertex, backbone_major, backbone_minor = make_SDKvertex(AA_molecule, atom_type1, atom_type2, AAbackbone, AAroots, AAbranches, AArings)

            print('backbone',backbone_vertex, 'backbone major', backbone_major, 'backbone minor', backbone_minor)
            print('ring', similar_ring, 'ring major',  major_ring, 'ring minor', minor_ring)
            major = backbone_major + major_ring
            minor = backbone_minor + minor_ring
            vertices = similar_ring + backbone_vertex

    else:
        similar_back, major_back, minor_back = make_SDKvertex(AA_molecule, atom_type1, atom_type2, AAbackbone, AAroots, AAbranches, AArings)
        major = major_back
        minor = minor_back
        vertices = similar_back

print('vertices:',vertices)
print('major:',major)
print('minor:',minor)

AA_contracted = contract_minor_nodes(minor, AA_molecule)
cg = list(AA_contracted.nodes())
# print(cg)
# print(AA_molecule.nodes())

sorted_node = []

CG_molecule = nx.Graph()

beads = make_bead(vertices, AA_nodes, atom_type1, atom_type2, AAbackbone, AA_molecule, major)
CG_nodes = update_nodes_aftergraining(vertices, AA_contracted, beads, AA_nodes)
print('Coarse Grained nodes',CG_nodes)

# for i in range(len(CG_nodes)):
#     print(i, CG_nodes[i].type)

for nd1,nd2 in enumerate(CG_nodes):
    cg_node1 = nd1
    cg_node2 = cg[(nd1 + 1) % len(cg)]
    if cg_node1 not in sorted_node and cg_node2 not in sorted_node:
        if AA_molecule.has_edge(cg[cg_node1], cg_node2)==True:
           continue
        elif len(AArings)!=0:
            for ring in AArings:
                if cg[cg_node1] in ring and cg_node2 in ring:
                    AA_contracted.add_edge(cg[cg_node1], cg_node2, weight=1)
                else:
                    for s,m in enumerate(major):
                        ts = m
                        ns = major[(s + 1) % len(major)]
                        if AA_contracted.has_edge(ts, ns) == True or ns==ts:
                            continue
                        elif ts in major_ring or ns in major_ring:
                            continue
                        else:
                            AA_contracted.add_edge(ts, ns, weight=1)
        elif len(AAbranches)!=0:
            for nd in cg:
                for c in cg:
                    if c != nd and AA_molecule.has_edge(nd, c - 1):
                        AA_contracted.add_edge(nd, c, weight=1)
        else:
            for s, m in enumerate(major):
                ts = m
                ns = major[(s + 1) % len(major)]
                if ns == ts:
                    continue
                else:
                    if AA_contracted.has_edge(ts, ns) == True:
                        continue
                    else:
                        AA_contracted.add_edge(ts, ns, weight=1)
            if len(AAbackbone) == 4:
                pass
            else:
                if AA_contracted.has_edge(major[0], major[-1]) == True and len(major) > 2:
                    AA_contracted.remove_edge(major[0], major[-1])
# print('AA back',AA_backbone)
CG_1 = nx.convert_node_labels_to_integers(AA_contracted)
draw_graph(AA_contracted)
draw_graph(CG_1)

print('CG_1')
for item1 in CG_1:
    print(item1, CG_nodes[item1], CG_nodes[item1].type, CG_nodes[item1].content)


###########################ITERATIONS#####################################
#First iteration
# CG_ring, CG_branches, CG_roots, CG_backbone = check_features(CG_1)






CG1_backbone = []
CG1_branches = []
test_smth = []

for item in AA_contracted.nodes:
    for AA_bak in AAbackbone:
        if item in AA_bak:
            CG1_backbone.append(item)
        else:
            CG1_branches.append(item)
            for n1 in CG1_branches:
                for n2 in CG1_branches:
                    if n1!=n2 and AA_contracted.has_edge(n1,n2) and sorted([n1,n2]) not in test_smth:
                        test_smth.append(sorted([n1,n2]))
                    else:
                        continue

print('Coarse Grain backbone', CG1_backbone)
print('Coarse Grain branches', test_smth)


# for bb in CG_backbone:
#     for i in bb:
#         print(i, CG_nodes[i], CG_nodes[i].type)
#
#
# if len(CG_branches)!=0:
#     if check_branch(CG_branches)== True:
#         CG_branch_vertex, CG_major_branch, CG_minor_branch = SDK_simple(CG_branches, CG_1)
#         print('First it vertices', CG_branch_vertex)
#         print('First it major', CG_major_branch)
#         print('First it minor', CG_minor_branch)
#     else:
#         print('check again')
#
#
# CG_1_contracted = contract_minor_nodes(CG_minor_branch, CG_1)
# print(CG_1_contracted.nodes())
# draw_graph(CG_1_contracted)
#
#
#
#
# CG_beads1 = make_bead(CG_branch_vertex,CG_nodes, atom_type1, atom_type2, CG_backbone, CG_1, CG_major_branch)
# print(CG_beads1)
# CG1_nodes = update_nodes_aftergraining(CG_branch_vertex, CG_1, CG_beads1, CG_nodes)
# print('First iteration', CG1_nodes)
#
# CG1_nodes = update_nodes_aftergraining(vertices, CG_1, beads, AA_nodes)
# print('First iteration', CG1_nodes)

# print(branch_vertex)
# print(branch_major)
# print(branch_minor)

#
# print('CG minor', CGminor)
# print('CGback_vertex', CGvertices)
# print('CG major', CGmajor)


# test_CG = nx.Graph()
#
# CGnodes_range = range(len(CG_nodes))
# print(CGnodes_range)
# CG_nds = []
# for i1, i2 in enumerate(CGnodes_range):
#     nd_1 = i1
#     bd_1 = CG_nodes[i1]
#     nd_2 = i2+1
#     bd_2 = CG_nodes[(i1+1)%len(CG_nodes)]
#     CG_nds.append(nd_1)
#     CG_nds.append(nd_2)

#
# CG_nds = list(dict.fromkeys(CG_nds))
# print(CG_nds)
# draw_graph(test_CG, CG_nds)



#######################################################################################
# contr = []
# for i in AA_contracted.nodes:
#     contr.append(i)
# print(contr)
# # print(check_features(AA_contracted))
#
# str = contr[0]
# stp = contr[-1]
# print(str, stp)
#
# CG_back = []
# CG_ring = nx.cycle_basis(AA_contracted)
# CG_branches = []
# print('AA_back',AA_backbone)
# for c in contr:
#     for AAbak in AA_backbone:
#         if c in AAbak:
#             CG_back.append(c)
#         elif len(CG_ring)==0 and c not in AAbak:
#              CG_branches.append([c])
#
# CG_backbone = []
# CG_backbone.append(CG_back)
# print('new back',CG_backbone)
# print('new ring',CG_ring)
# print('new branch',CG_branches)
# #
# draw_graph(AA_contracted)


# print(SDK_small(AA_contracted, roots))
# print(make_SDKvertex(CG_molecule, atom_type1, atom_type2, CG_backbone,roots, CG_branches, CG_ring))
#
# print(CG_branch)
# print(nx.dfs_tree(AA_contracted))
# print(make_SDKvertex(AA_contracted, atom_type1, atom_type2, AA_backbone, roots, sad, rings))
# print(AA_contr)
# print(type(AA_contr))
# print(AA_contr[0], len(AA_contr)-1)
#
# start = list(AA_contracted.nodes())[0]
# stop = len(AA_contracted) - 1

#iterations


#####output
# make_xyz('CG_branches.xyz',CG_nodes)
# make_xyz('AA_branches.xyz', AA_nodes)

# draw_graph(AA_molecule)
