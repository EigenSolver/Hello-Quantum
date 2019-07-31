# -*- coding: utf-8 -*-
# File name: problem2.py

# Huawei Quantum Programming Competition
# Author: QuantumSupremacy Group


import numpy as np
from scipy.optimize import minimize
from projectq import MainEngine
from projectq.ops import H, Rx, Ry, Rz, Measure, All, CNOT
from projectq.cengines import (MainEngine,
                               AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               InstructionFilter,
                               DecompositionRuleSet)
import projectq.setups.decompositions as rules

from projectq.backends import Simulator
from hiq.projectq.backends import SimulatorMPI
from hiq.projectq.cengines import GreedyScheduler, HiQMainEngine

def bitstring(n):
    if n > 1:
        return ['0'+bits for bits in bitstring(n-1)]+['1'+bits for bits in bitstring(n-1)]
    elif n == 1:
        return ['0', '1']
    else:
        print('error')

def parse_mapper(mapper):
    """
    Args:
        mapper(string):
    Returns:
        n(int): number of qubits
        edges(list[tuple]): list of 2 elements tuples representing the edges in the mapper, indices start from zero
    """
    info=mapper.split("\n")
    return int(info[0]), list(map(lambda x: tuple(map(lambda y: int(y)-1, x.split(','))),info[1:]))

def reduce_connection(n,edges):
    """
    cut the extra edges and reduce the coupling between qubits

    Args:
        n(int): number of qubits
        edges(list[tuple]): list of 2 elements tuples representing the edges in the graph
    Returns:
        reduced_edges(list[tuple]): list of 2 elements tuples representing the edges in the tree
    """

    reduced_edges=edges.copy()
    degree_list=np.zeros(n)
    for edge in edges:
        degree_list[edge[0]]+=1
        degree_list[edge[1]]+=1

    for edge in edges:
        if degree_list[edge[0]]>2 and degree_list[edge[1]]>2:
            reduced_edges.remove(edge)
            degree_list[edge[0]]-=1
            degree_list[edge[1]]-=1

    return reduced_edges

def cost_func(param, target_state, p, eng, connections=None):
    '''
    VQE core part
    Args:
        target_state(string): bitstring contains of 0/1
        p(int): circuit depth of VQC
        eng(MainEngine): engine of simulator
    Return:
        cost(int): cost function 
    '''

    tolerance=0.001
    assert np.abs(sum([np.abs(i)**2 for i in target_state])-1)<tolerance, "Invalid Target State"

    n=int(np.log2(len(target_state))) 
    assert len(param)==3*n*p
    param=param.reshape((n,3*p))

    count = 0
    all_string = bitstring(n)


    qureg = eng.allocate_qureg(n)
    All(H) | qureg
    for m in range(p):
        for i in range(n):
            Rx(param[i, 3*m+0]) | qureg[i]
            Ry(param[i, 3*m+1]) | qureg[i]
            Rz(param[i, 3*m+2]) | qureg[i]
        
        if connections is None:
            CNOT | (qureg[n-1], qureg[0])
            for i in range(n-1):
                CNOT | (qureg[i], qureg[i+1])
        else:
            for edge in connections:
                CNOT | (qureg[edge[0]], qureg[edge[1]])
                
    eng.flush()
    
    for i in range(len(target_state)):
        count += target_state[i]*eng.backend.get_amplitude(all_string[i], qureg)
    All(Measure) | qureg

    return 1-np.abs(count)


def circuit_generator(eng, target_state, mapper):
    """
    Args:
        eng:
        The engine type should be defined in "__main__" function.
        target_state(list):
        A target quantum state vector like: [0.5. 0.5. 0.5. 0.5]. And this
        simple example can be prepared from 'H | qubit[0]; H | qubit[1]'
        where the 'qubit' are in a |0> state at very first.
        mapper(string):
        The mapper of a real quantum chip like: '3\n1,2\n2,3'. This simple
        example of mapper means that there are 3 qubits in this given chip
        and you can play 'CNOT' on 'qubit1' and 'qubit2'. You can also play
        'CNOT' on 'qubit2' and 'qubit3'. But you can not play 'CNOT' on
        'qubit1' and 'qubit3' because '1,3' is not in the given mapper.
    Returns:
        simulated_circuit(string):
        After your calculation you may get a final quantum circuit like:
        'H | qubit[0]; H | qubit[1]'. Each step seperated by a '; '. In
        our score we will search operations and qubit index in your result
        and run it on our backend. We will compare the final state your
        circuit produces with our target state by ourselves.
    """
    n=int(np.log2(len(target_state))) 
    p=1 
    param = np.ones((n, 3*p)).flatten()
    assert n==int(mapper[0])
    connections=reduce_connection(*parse_mapper(mapper))
    
    result=minimize(cost_func,x0=param, args=(target_state,p,eng,connections),method='Nelder-Mead')
    opt_param=result.x.reshape((n,3*p)) # get optimal parameter
    
    simulated_circuit = ''
    for i in range(n):
        simulated_circuit+='H | qubit[{}]; '.format(i)
    for m in range(p):
        for i in range(n):
            simulated_circuit+="Rx({0}) | qubit[{1}]; ".format(opt_param[i, 3*m+0],i)
            simulated_circuit+="Ry({0}) | qubit[{1}]; ".format(opt_param[i, 3*m+1],i)
            simulated_circuit+="Rz({0}) | qubit[{1}]; ".format(opt_param[i, 3*m+2],i)
        for edge in connections:
                simulated_circuit+="CNOT | (qubit[{0}], qubit[{1}]); ".format(*edge)
        
    return simulated_circuit[:-2]

if __name__=="__main__": # need to set


    backend=Simulator(gate_fusion=True)
    cache_depth = 10
    rule_set = DecompositionRuleSet(modules=[rules])
    engines = [TagRemover(), LocalOptimizer(cache_depth), AutoReplacer(rule_set), TagRemover(), LocalOptimizer(cache_depth)
               ]
    
    eng= MainEngine(backend, engines)
    
    target_state = [0.5, 0.5, 0.5, 0.5]
    mapper = '2\n1,2'
#    
#    backend = SimulatorMPI(gate_fusion=True, num_local_qubits=int(np.log2(len(target_state))))
#
#    cache_depth = 10
#    rule_set = DecompositionRuleSet(modules=[rules])
#    engines = [TagRemover(), LocalOptimizer(cache_depth), AutoReplacer(rule_set), TagRemover(), LocalOptimizer(cache_depth), GreedyScheduler()
#            ]
#    eng = HiQMainEngine(backend, engines)
    

    circuit = circuit_generator(eng, target_state, mapper)

    print(circuit)
