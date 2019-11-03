import networkx as nx
import matplotlib.pyplot as plt

def draw_graph(graph, nodes=None, others=None):
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