import networkx as nx
import random
import matplotlib.pyplot as plt

def graph(num_qubits):
    
    """
    Generate a random unidirectional graph with automated source and destination nodes and colors.
    
    Args:
        num_qubits (int): Number of nodes in the graph.

    Returns:
        nx.DiGraph: The generated directed graph.
        list of tuples: List of edges as two-tuple list.
        list of floats: List of corresponding weights for each edge.
        int: Source node.
        int: Destination node.
        nodes: no. of nodes
    """
    # Calculate the number of edges based on the density factor
    num_edges = random.randint(num_qubits, (num_qubits * (num_qubits - 1))/2)

    # Create an empty directed graph (DiGraph)
    G = nx.DiGraph()

    # Add nodes to the graph
    G.add_nodes_from(range(num_qubits))

    # Initialize lists to store edges and weights
    edges_list = []
    weights_list = []

    # Create a list of all possible edge pairs
    edge_pairs = [(i, j) for i in range(num_qubits) for j in range(num_qubits) if i != j]

    # Shuffle the list of edge pairs
    random.shuffle(edge_pairs)

    # Add random weighted directed edges to the graph while avoiding bidirectional edges and self-loops
    for i in range(num_edges):
        source_node, destination_node = edge_pairs[i]

        # Avoid bidirectional edges and self-loops
        if not G.has_edge(destination_node, source_node):
            weight = random.uniform(0,1)

            # Add the edge to the graph
            G.add_edge(source_node, destination_node, weight=weight)

            # Append the edge as a two-tuple and its weight to the respective lists
            edges_list.append((source_node, destination_node))
            weights_list.append(weight)

    # Assign automated source and destination nodes
    source_node = random.randint(0, num_qubits - 1)
    destination_node = random.randint(0, num_qubits - 1)
    
    # Ensure source and destination nodes are different
    while source_node == destination_node:
        destination_node = random.randint(0, num_qubits - 1)

    if not nx.has_path(G, source_node, destination_node):
        # If there is no path, create an edge from source to destination
        weight = random.uniform(0,1)
        G.add_edge(source_node, destination_node, weight=weight)
        edges_list.append((source_node, destination_node))
        weights_list.append(weight)

    nodes = G.number_of_nodes()

    return G, edges_list, weights_list, source_node, destination_node, nodes

# Graph Details:
num_qubits = 8

G, edges_list, weights_list, source_node, destination_node, nodes = graph(num_qubits)

# Create a color mapping for nodes
node_colors = ['green' if node == source_node else 'blue' if node == destination_node else 'skyblue' for node in G.nodes()]

# Visualize the graph with node colors (optional)
layout = nx.spring_layout(G)
nx.draw(G, pos=layout, with_labels=True, node_size=300, node_color=node_colors, font_size=10, font_color='black', font_weight='bold')
edge_labels = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G, pos=layout, edge_labels=edge_labels)
plt.show()

# Print the generated lists and nodes
print("List of Edges (as two-tuple list):")
print(edges_list)
print("List of Weights:")
print(weights_list)
print("Source Node:", source_node)
print("Destination Node:", destination_node)