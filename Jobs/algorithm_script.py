#!/usr/bin/env python
# coding: utf-8

def Static_Model(edges, weights, nodes, source, destination):

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
    import string
    import time

    # AWS imports: Import Braket SDK modules
    from braket.circuits import Circuit, Gate, Observable
    from braket.devices import LocalSimulator
    from braket.aws import AwsDevice, AwsQuantumTask

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

            quadratic_coeff = quadratic.get(('x['+str(qubit_1)+']', 'x['+str(qubit_2)+']'))

            circuit.zz(qubit_1, qubit_2, quadratic_coeff*gamma)

        ## Implementing the Mixer Hamiltonian ##
        for qubit in range(n_qubits):
            circuit.rx(qubit, 2*beta) # theta=2*beta because rx rotates the qubit about X by theta/2 angle

        return circuit


    # * <b>Remember that there are 3 parameters that can be varied - beta, gamma and penalty p </b>

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

        expectation /= 1000
        return expectation

    ## Now we measure the circuit ##
    def expectation_execute_circuit(param):
        ## Set up the device to run the circuit
        device = LocalSimulator()

        ## QAOA parameters to be optimized such that the eigenvalue Cost(β, γ) can be minimized ##
        beta = param[0]
        gamma = param[1]

        circuit = create_circuit(beta, gamma)

        shots = 1000
        result = device.run(circuit, shots).result()
        counts = result.measurement_counts

        return compute_expectation(counts, shots)    

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
    result = device.run(circuit, shots = 1000).result()
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

    i = 0
    max = 0
    for i in range(len(possible_states)):
        if counts[possible_states[i]]>counts[possible_states[max]]:
            max = i
        #print(possible_states[i]+':'+str(counts[possible_states[i]]))

    return possible_states, counts

def start_here():
    edges = [(1,2),(1,3),(2,3),(2,4),(2,5),(3,4),(3,5),(4,5),(4,6),(5,6)]
    weights = [5,8,2,7,4,3,1,1,2,10]

    nodes = 6
    ## Set the source and destination nodes ##
    source = 1
    destination = 6
    possible_states, counts = Static_Model(edges, weights, nodes, source, destination)
    
    print(possible_states)
    print(counts)

start_here()