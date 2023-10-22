#!/usr/bin/env python
# coding: utf-8

# # Graph #

# In[1]:


def graph(num_qubits):
    import networkx as nx
    import random
    import matplotlib.pyplot as plt
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

# Example usage:
'''
# Create a color mapping for nodes
node_colors = ['green' if node == source else 'blue' if node == destination else 'skyblue' for node in G.nodes()]

# Visualize the graph with node colors (optional)
layout = nx.spring_layout(G)
nx.draw(G, pos=layout, with_labels=True, node_size=300, node_color=node_colors, font_size=10, font_color='black', font_weight='bold')
edge_labels = {(u, v): d['weight'] for u, v, d in G.edges(data=True)}
nx.draw_networkx_edge_labels(G, pos=layout, edge_labels=edge_labels)
plt.show()

# Print the generated lists and nodes
print("List of Edges (as two-tuple list):")
print(edges)
print("List of Weights:")
print(weights)
print("Source Node:", source)
print("Destination Node:", destination)
'''


# # QAOA Algorithm #

# In[2]:


def Static_Model(G, edges, weights, source, destination, nodes):
    from pyqubo import Array,Constraint,Placeholder
    
    x = Array.create(name = 'x', shape =len(edges), vartype = 'BINARY') 

    i = 0
    fcost = 0
    for i in range(len(edges)):
        fcost += (weights[i]*x[i]) ## Put the expression for constraints in

    p = 27
    node=1
    while(node<=nodes):
        sum = 0
        for j in range(len(edges)):
            if edges[j][0] == node: ##Outgoing from the current node
                sum += x[j]
            elif edges[j][1] == node: ##Incoming into the current node
                sum -= x[j]

        if node == source:
            sum -= 1
        elif node== destination:
            sum += 1

        fcost += p*((sum)**2)
        node+=1

    model = fcost.compile()

    linear, quadratic, offset = model.to_ising()

    # general imports
    import numpy as np
    import matplotlib.pyplot as plt
    # magic word for producing visualizations in notebook
    #get_ipython().run_line_magic('matplotlib', 'inline')
    import string
    import time

    # AWS imports: Import Braket SDK modules
    from braket.circuits import Circuit, Gate, Observable
    from braket.devices import LocalSimulator
    from braket.aws import AwsDevice, AwsQuantumTask


    # In[28]:


    def create_circuit(beta, gamma):
        ## initializing the initial qubit state with H gates ##
        circuit = Circuit()
        n_qubits = len(edges)

        for qubit in range(n_qubits):
            circuit.h(qubit)

        ## Implementing the problem Hamiltonian ##
        for qubit in range(n_qubits):
            linear_coeff = linear.get('x['+str(qubit)+']')
            circuit = circuit.rz(qubit, -1*linear_coeff)

        #Algorithmic method to add the ZZ gates - CHECK TO SEE IF IT AFFECTS THE RESULTS(it should'nt because they commute)
        for i in range(len(quadratic)):
            qubit_1 = int(list(quadratic.keys())[i][0][2])
            qubit_2 = int(list(quadratic.keys())[i][1][2])
            key = ('x[' + str(qubit_1) + ']', 'x[' + str(qubit_2) + ']')
        
            # Check if the key exists in the dictionary
            if key in quadratic:
                quadratic_coeff = quadratic[key]
                # The Ising-Coupling Gate
                #circuit.zz(qubit_1, qubit_2, quadratic_coeff*gamma)
                circuit.cnot(control=qubit_1, target=qubit_2)
                circuit.rz(qubit_2, quadratic_coeff*gamma)
                circuit.cnot(control=qubit_1, target=qubit_2)
                
        ## Implementing the Mixer Hamiltonian ##
        for qubit in range(n_qubits):
            circuit.rx(qubit, 2*beta) # theta=2*beta because rx rotates the qubit about X by theta/2 angle

        return circuit


    # * <b>Remember that there are 3 parameters that can be varied - beta, gamma and penalty p </b>

    # In[29]:


    ## Expectation value of the Hamiltonian is basically the expected cost value which we can get from an average of the
    ## cost values over all states that have occurred ##
    def compute_expectation(counts, shots):

        expectation = 0
        sum = 0
        states = list(counts.keys())
        for i in range(len(states)):
            state = states[i] # string variable of the current qubit states
            state_cost = 0
            for j in range(len(state)): # Convention of the states is that the left-most qubit is the first qubit q0
                state_cost = state_cost + int(state[j])*weights[j]

            expectation = expectation + state_cost*counts.get(state)

        expectation /= 10000
        return expectation


    # In[30]:


    ## Now we measure the circuit ##
    def expectation_execute_circuit(param):
        ## Set up the device to run the circuit
        device = LocalSimulator()

        ## QAOA parameters to be optimized such that the eigenvalue Cost(β, γ) can be minimized ##
        beta = param[0]
        gamma = param[1]

        circuit = create_circuit(beta, gamma)

        shots = 10000
        result = device.run(circuit, shots).result()
        counts = result.measurement_counts

        return compute_expectation(counts, shots)    


    # In[31]:


    ## Classical Optimizer ##
    from scipy.optimize import minimize

    res = minimize(expectation_execute_circuit,
                   [1.0, 1.0],
                   method='COBYLA')

    # ## Analyzing the Results
    beta = res.get('x')[0]
    gamma = res.get('x')[1]
    circuit = create_circuit(beta, gamma)

    device = LocalSimulator()
    result = device.run(circuit, shots = 10000).result()
    counts = result.measurement_counts

    # plot using Counter
    #plt.bar(counts.keys(), counts.values())
    #plt.xlabel('bitstrings')
    #plt.ylabel('counts')


    # ## Post-Processing
    # 
    # Remove the output states that are not possible (among the top ten most probable states) and then check for the most probable states.
    # Check for joined paths, that is check that if a path enum ending with a number exists, then another path enum starting with the same number also exists, unless its a source or destination. They should always be 1.
    def check_state(s):

        ## Firstly check if the path starts from a source and ends at a destination
        source_flag = False
        destination_flag = False
        multiple_branches = False
        continuity_flag = True

        starting_nodes = []
        ending_nodes = []

        ## Check to see if the source and destination nodes exist in the network
        i=0
        for i in range(len(s)):
            if(s[i] == '1'):

                if(edges[i][0] == source):
                    source_flag = True
                if(edges[i][1] == destination):
                    destination_flag = True

                starting_nodes += [edges[i][0]]
                ending_nodes += [edges[i][1]]

        ## Now check if a node repeats itself in starting or ending_nodes. If yes, set multiple_branches
        i = 0
        for i in range(len(starting_nodes)):
            cnt1 = starting_nodes.count(starting_nodes[i])
            cnt2 = ending_nodes.count(ending_nodes[i])
            if cnt1 > 1 or cnt2 > 1:
                multiple_branches = True
                break
        ## Then iteratively go through ending nodes and check if the same node exists in the next starting_nodes index                
        ## This is an easier check for continuity and necessarily requires the edges nodes to be in some order
        ## Also go with the thumb rule that the destination node will be the last value of ending_nodes
        for i in range(len(ending_nodes)-1):
            if starting_nodes[i+1] != ending_nodes[i]:
                continuity_flag = False
                break
        if source_flag and destination_flag and continuity_flag and (not multiple_branches):
            return True
        else:
            return False

    states = list(counts.keys())
    possible_states = []
    i = 0
    for i in range(len(states)):
        s = states[i]
        flag = check_state(s)
        if flag:
            possible_states += [s]
    del states
    
    # The cost of the optimized Path:
    count_states = []
    for state in possible_states:
        if state in counts:
            count_states.append(str(counts[state]))
        else:
            count_states.append('0')  # Use 0 if the state is not in counts
    
    dict_states = dict(zip(possible_states, count_states))
    
    if dict_states:
        max_state = max(dict_states, key=lambda x: int(dict_states[x]))
        cost_opt = np.sum(np.multiply([int(char) for char in max_state], weights))
    else:
        cost_opt = 0  # Handle the case where dict_states is empty
    
    return cost_opt


# # Approximation Ratio #

# In[3]:


#Approximation Ratio:
def approx_ratio(G, edges, weights, source, destination, nodes):
    import networkx as nx
    # Find the most efficient path using Dijkstra's algorithm
    actual_cost_opt = nx.dijkstra_path_length(G, source=source, target=destination)
    cost_opt = Static_Model(G, edges, weights, source, destination, nodes)
    if cost_opt == 0:
        actual_cost_opt, aprox_ratio = approx_ratio(G, edges, weights, source, destination, nodes)
    else:
        aprox_ratio = actual_cost_opt/cost_opt
    return [actual_cost_opt,aprox_ratio]

# Average Approximation Ratio:
def av_approx_ratio(no_of_iters,G, edges, weights, source, destination, nodes):
    sum_approx_ratio = 0
    for i in range(no_of_iters):
        '''
        if approx_ratio != 0:
            approx_ratio_value = approx_ratio(G, edges, weights, source, destination, nodes)[1]
            sum_approx_ratio += approx_ratio_value
        else:
            i = i - 1
        '''
        approx_ratio_value = approx_ratio(G, edges, weights, source, destination, nodes)[1]
        sum_approx_ratio += approx_ratio_value
    av_approx_ratio = sum_approx_ratio/no_of_iters
    return av_approx_ratio

def test_1(no_of_iters,G, edges, weights, source, destination, nodes):
    import matplotlib.pyplot as plt
    iter_array = list(range(1,no_of_iters+1))
    av_approx = []
    approx = []
    actual = []
    for i in iter_array:
        '''
        if approx_ratio(G, edges, weights, source, destination, nodes)[1] != 0:
            av_approx.append(av_approx_ratio(i,G, edges, weights, source, destination, nodes))
            approx.append(approx_ratio(G, edges, weights, source, destination, nodes)[1])
            #actual.append(approx_ratio(G, edges, weights, source, destination, nodes)[0])
        else:
            i = i - 1
        '''
        av_approx.append(av_approx_ratio(i,G, edges, weights, source, destination, nodes))
        approx.append(approx_ratio(G, edges, weights, source, destination, nodes)[1])
        actual.append(approx_ratio(G, edges, weights, source, destination, nodes)[0])
    plt.plot(iter_array, av_approx,'o-')
    plt.plot(iter_array, approx,'o--')
    plt.plot(iter_array, actual)
    plt.show()

'''
def test_2(num_of_qubits,G, edges, weights, source, destination, nodes):
    import matplotlib.pyplot as plt
    qubit_array = list(range(3,no_of_qubits+1))
    av_approx = []
    approx = []
    actual = []
    for i in qubit_array:
        if approx_ratio(G, edges, weights, source, destination, nodes)[1] != 0:
            av_approx.append(av_approx_ratio(i,G, edges, weights, source, destination, nodes))
            approx.append(approx_ratio(G, edges, weights, source, destination, nodes)[1])
            #actual.append(approx_ratio(G, edges, weights, source, destination, nodes)[0])
        else:
            i = i - 1
        av_approx.append(av_approx_ratio(i,G, edges, weights, source, destination, nodes))
        approx.append(approx_ratio(G, edges, weights, source, destination, nodes)[1])
        actual.append(approx_ratio(G, edges, weights, source, destination, nodes)[0])
    plt.plot(iter_array, av_approx,'o-')
    plt.plot(iter_array, approx,'o--')
    plt.plot(iter_array, actual)
    plt.show()
'''


# # Testing of the Model #

# In[4]:


def test_model(num_qubits,no_of_iters):
    G, edges, weights, source, destination, nodes = graph(num_qubits)
    test_1(no_of_iters,G, edges, weights, source, destination, nodes)


# In[6]:


for i in range(5,7):
    num_qubits = i
    no_of_iters = 10
    test_model(num_qubits,no_of_iters)


# In[ ]:




