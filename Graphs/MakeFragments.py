import numpy as np
import numpy.linalg
from rdkit import Chem
import rdkit
from rdkit.Chem import AllChem
import networkx as nx

class Bead:
#     'Collection of double bonded beads'
    def __init__(self, bead_nodes, bead_CoM,  bead_mass, bead_id):
        self.content = bead_nodes
        self.coords = bead_CoM
        self.weight = bead_mass
        self.type = bead_id


def simple_fragment(nodes):
    vertices = []
    for c in range(0, len(nodes), 3):
        chn = nodes[c: c + 3]
        vertices.append(chn)
    return vertices

def heavy_fragment(fragnodes, othernodes, graph):
    edge_frag = []
    back_frag = []
    for num, at in enumerate(fragnodes):
        at1 = at
        at2 = fragnodes[(num + 1) % len(fragnodes)]
        if at2 == fragnodes[0] and at1 == fragnodes[-1]:
            continue
        else:
            back_frag.append(nx.shortest_path(graph, source=at1, target=at2))

    if graph.has_edge(fragnodes[0], othernodes[0]):
        edge_frag.append([fragnodes[0], othernodes[0]])
    if graph.has_edge(fragnodes[-1], othernodes[-1]):
        edge_frag.append([fragnodes[-1], othernodes[-1]])

    groups = sorted(back_frag + edge_frag)

    vertices = []
    smth = []
    for v1 in groups:
        for v2 in groups:
            v1 = set(v1)
            v2=set(v2)
            if v1!=v2 and v1.issubset(v2):

                smth.append(list(v1))
    print(smth)
    for group in groups:
        if sorted(group) not in smth:
            vertices.append(sorted(group))

    return back_frag, edge_frag, vertices

def rigid_fragment(graph, in_list):
    vertices = []
    print(in_list)
    Lap = nx.normalized_laplacian_matrix(graph, nodelist=None)
    for u, v in enumerate(in_list):
        n1 = v
        n2 = in_list[(u+1)%len(in_list)]
        sim_rank = abs(Lap.A[n1][n2])
        if sim_rank == 0 or sim_rank == 1:
            continue
        elif sim_rank > 0.5 and sorted([n1, n2]) not in vertices:
            vertices.append([n1, n2])
    print(vertices)
    return vertices



def branch_fragment(backbone, roots, graph):
    short_back = []
    for atom in backbone:
        if atom in roots:
            continue
        else:
            short_back.append(atom)
    back_frag, edge_frag, vertices = heavy_fragment(short_back, roots, graph)

    return vertices

def linear_fragment(backbone, others, graph):
    vertices = []
    groups1 = []

    if len(others)==0:
        if len(backbone)<=3:
            vertices.append(backbone)
        elif len(backbone) == 4:
            for c in range(0, len(backbone), 2):
                chn = backbone[c: c + 2]
                vertices.append(chn)
        elif len(backbone)>4:
            sb = []
            if len(backbone) % 3 == 0:
                pass
            else:
                sb.append([backbone[0], backbone[1]])
                for s in sb:
                    groups1.append(s)
                    for b in s:
                        backbone.remove(b)
            groups2 = simple_fragment(backbone)
            vertices = groups1 + groups2
    else:
        back_frag, edge_frag, vertices = heavy_fragment(others, backbone, graph)
        print(vertices)

        print(others, backbone)
    return vertices

def rank_nodes(vertex, centrality, graph):
    major = []
    minor = []
    touched = []

    cluster = list(vertex.nodes())
    # print(cluster)

    for i, nodes in enumerate(cluster):
        nd1 = nodes
        nd2 = cluster[(i + 1) % len(cluster)]
        if round(centrality[nd1], 3) == round(centrality[nd2], 3):
            touched.append(nd1)
            touched.append(nd2)
        elif round(centrality[nd1], 3) > round(centrality[nd2], 3):
            major.append(nd1)
            minor.append(nd2)
        elif round(centrality[nd1], 3) < round(centrality[nd2], 3):
            major.append(nd2)
            minor.append(nd1)

        if nd1 not in major and nd1 not in minor and nd2 not in major and nd2 not in minor and nd1!=nd2:
            major.append(nd1)
            minor.append(nd2)

        for m in major:
            for j in major:
                if m!=j and m in cluster and j in cluster:
                    major.remove(m)
                    minor.append(m)

    major = list((dict.fromkeys(major)))
    minor = list((dict.fromkeys(minor)))

    return major, minor


def make_fragment(mol_graph, mol_nodes, rings, branches, roots, backbone, rigid_ring, rigid_branch, rigid_backbone, others):
    vertices = []
    nodes_list = list(mol_graph.nodes())
    if len(rings)==0 and len(branches)==0 and rigid_backbone == False:
        vertices = linear_fragment(backbone, others, mol_graph)
        # print(vertices)
    elif len(rings)==0 and len(branches)==0 and rigid_backbone == True:
        vertices = rigid_fragment(mol_graph, backbone)
        print(vertices)
    elif len(rings)!=0:
        if backbone == others or len(backbone)==0:
            vertices = rings
            print(vertices)
        elif rigid_ring == True:
            vertices_ring = rings
            vertices_backbone = branch_fragment(backbone, roots, mol_graph)
            vertices = vertices_ring + vertices_backbone
            print(vertices)
    elif len(branches)!=0:
        if rigid_branch == False:
            for branch in branches:
                if len(branch)==1:
                    vertices = simple_fragment(nodes_list)
    major_nd=[]
    minor_nd=[]
    for vertex in vertices:
        sub_graph = nx.Graph(mol_graph.subgraph(vertex))
        sub_centrality = nx.pagerank_numpy(sub_graph, alpha=0.85, weight='w')
        mjr_nd, mnr_nd = rank_nodes(sub_graph, sub_centrality, mol_graph)
        for m in mjr_nd:
            major_nd.append(m)
        for n in mnr_nd:
            minor_nd.append(n)

    return vertices, major_nd, minor_nd

def contract(minor, major, graph):
    FragGraph = nx.Graph(graph)
    for nd in minor:
        if nd in FragGraph.nodes():
            FragGraph.remove_node(nd)
        else:
            continue
    all_nodes = list(graph.nodes())
    FragNodes = FragGraph.nodes()
    # print(sorted(major))
    for m, mj in enumerate(major):
        ts = mj
        ns = major[(m + 1) % len(sorted(major))]
        print(ts, ns)
        if FragGraph.has_edge(ts, ns) == True:
            continue
        elif graph.has_edge(all_nodes[0], ts-1) and all_nodes[0] in FragNodes:
           FragGraph.add_edge(all_nodes[0], ts, weight=1)

        FragGraph.add_edge(ts, ns, weight=1)
        if FragGraph.has_edge(major[0], major[-1]) == True and len(major) > 2:
            FragGraph.remove_edge(major[0], major[-1])

    return FragGraph, FragNodes









