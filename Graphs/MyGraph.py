import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg

class MyGraph():
    ''' a wrapper for a networkx graph, and related rings/branches/roots/backbone to perform a single iteration of contraction '''
    def __init__(self, graph, others=[]):
        ''' perform initial analysis 
        :graph: input should be a networkx.graph object '''
        self.graph = graph
        self.others = others
         # initialize rings/branches/roots/backbone
        rings, branches, roots, backbone = self._find_regions(self.graph)
        rigid_ring, rigid_branch, rigid_backbone = self._check_rigid(self.graph, rings, branches, backbone)
        
        self.rings = rings 
        self.branches = branches
        self.roots = roots
        self.backbone = backbone
        self.rigid_ring = rigid_ring
        self.rigid_branch = rigid_branch
        self.rigid_backbone = rigid_backbone

    def _nodes2vertices(self, similarity):
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

    def _find_rings(self, graph, rings):
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

    def _find_branches(self, graph):
        graph_branches = []
        test_brnh = []
        temp_branches = []
        roots = []
        start = list(graph.nodes())[0]
        # stop = len(graph) - 1
        stop = max(graph.nodes)

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
                    graph_branches = self._nodes2vertices(test_brnh)
                else:
                    for root in roots:
                        if graph.has_edge(nd1, root) and graph.degree(nd1) == 1 and [nd1] not in graph_branches:
                            graph_branches.append([nd1])

        return graph_branches, roots, backbone

    def _find_regions(self, graph):

        matches_rings = nx.cycle_basis(graph)
        if len(matches_rings)!=0:
            rings, roots, backbone = self._find_rings(graph, matches_rings)
            branches = []
        else:
            branches, roots, backbone = self._find_branches(graph)
            rings = []

        return rings, branches, roots, backbone

    def _find_rigid(self, graph, region):
        conj_func = []

        for (u, v, d) in graph.edges(data=True):
            if u in region and v in region:
                if graph.has_edge(u, v):

                    if graph[u][v]['weight'] == 1:
                        conj_func.append(-1)
                    elif graph[u][v]['weight'] > 1:
                        conj_func.append(1)

        return conj_func

    def _conj_region(self, conj):

        conj_func = sum(conj)
        if conj_func>0:
            rigid=True
        else:
            rigid=False
        return rigid

    def _check_rigid(self, graph, rings, branches, backbone):
        rigid_ring = 0
        rigid_branch = 0
        rigid_backbone = 0

        if len(rings)!=0:
            for ring in rings:
                ring_conj_func = self._find_rigid(graph, ring)
            rigid_ring = self._conj_region(ring_conj_func)
        elif len(branches)!=0:
            for branch in branches:
                branch_conj_func = self._find_rigid(graph, branch)
            rigid_branch = self._conj_region(branch_conj_func)
        else:
            backbone_conj_func = self._find_rigid(graph, backbone)
            rigid_backbone = self._conj_region(backbone_conj_func)

        return rigid_ring, rigid_branch, rigid_backbone

    def make_fragment(self):
        vertices = []
        if len(self.rings)==0 and len(self.branches)==0 and self.rigid_backbone == False:
            vertices = self._linear_fragment(self.backbone, self.others, self.graph)
            print(vertices)
        elif len(self.rings)!=0:
            if self.backbone == self.others or len(self.backbone)==0:
                vertices = self.rings
                print(vertices)
            elif self.rigid_ring == True:
                vertices_ring = self.rings
                vertices_backbone = self._branch_fragment(self.backbone, self.roots, self.graph)
                vertices = vertices_ring + vertices_backbone
                print(vertices)
        elif len(self.branches)!=0:
            if self.rigid_branch == False:
                for branch in self.branches:
                    if len(branch)==1:
                        vertices = self.simple_fragment(self.mol_nodes)
        major_nd=[]
        minor_nd=[]
        for vertex in vertices:
            sub_graph = nx.Graph(self.graph.subgraph(vertex))
            sub_centrality = nx.pagerank_numpy(sub_graph, alpha=0.85, weight='w')
            mjr_nd, mnr_nd = self._rank_nodes(sub_graph, sub_centrality, self.graph)
            for m in mjr_nd:
                major_nd.append(m)
            for n in mnr_nd:
                minor_nd.append(n)

        return vertices, major_nd, minor_nd

    def _simple_fragment(self, nodes):
        vertices = []
        for c in range(0, len(nodes), 3):
            chn = nodes[c: c + 3]
            vertices.append(chn)
        return vertices

    def _heavy_fragment(self, fragnodes, othernodes, graph):
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

    def _branch_fragment(self, backbone, roots, graph):
        print(backbone, roots)
        short_back = []
        for atom in backbone:
            if atom in roots:
                continue
            else:
                short_back.append(atom)
        back_frag, edge_frag, vertices = self._heavy_fragment(short_back, roots, graph)
        return vertices

    def _linear_fragment(self, backbone, others, graph):
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
                groups2 = self._simple_fragment(backbone)
                vertices = groups1 + groups2
        else:
            back_frag, edge_frag, vertices = self._heavy_fragment(others, backbone, graph)
            print(vertices)

            print(others, backbone)
        return vertices

    def _rank_nodes(self, vertex, centrality, graph):
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

    def contract(self, minor, major):
        graph=self.graph
        FragGraph = nx.Graph(graph)
        for nd in minor:
            if nd in FragGraph.nodes():
                FragGraph.remove_node(nd)
            else:
                continue
        all_nodes = list(graph.nodes())
        FragNodes = FragGraph.nodes()
        print(sorted(major))
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

    def draw_graph(self, nodes=None, others=None):
        graph = self.graph
        color_map = []
        atom_type2 = []

        esingle = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 1]
        etriple = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 3]
        edouble = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 2]
        earom = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 1.5]
        # var = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] == 0.75]

        pos = nx.spring_layout(graph)
        nx.draw_networkx_nodes(graph, pos, node_size=500)
        # nx.draw_networkx_nodes(graph, pos, node_size=250, nodelist=nodes, node_color='g')

        if nodes == None:
            nodes=graph.nodes()
        else:
            for node in range(len(nodes)):
                # print(node, nodes[node].type)

                if nodes[node].type == 'O':

                    nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=others, node_color='r')
                elif nodes[node].type == 'Cl':

                    nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=others, node_color='y')
                elif nodes[node].type == 'F':

                    nx.draw_networkx_nodes(graph, pos, node_size=500, nodelist=others, node_color='g')

        nx.draw_networkx_edges(graph, pos, edgelist=etriple, width=12, alpha=0.5, edge_color='b', node_color=color_map)
        nx.draw_networkx_edges(graph, pos, edgelist=edouble, width=12, alpha=0.5, edge_color='g', node_color=color_map)
        nx.draw_networkx_edges(graph, pos, edgelist=esingle, width=10, node_color=color_map)
        nx.draw_networkx_edges(graph, pos, edgelist=earom, width=12, alpha=0.5, edge_color='r', node_color=color_map)
        # nx.draw_networkx_edges(graph, pos, edgelist=var, width=3, style='dashed')
        nx.draw_networkx_labels(graph, pos, font_size=20, font_family='sans-serif', node_color=color_map)
        plt.axis('off')
        plt.show()