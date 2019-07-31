# Huawei Quantum Programming Competition
# Author: QuantumSupremacy Group

import numpy as np

import projectq as pq
from projectq import MainEngine
from projectq.ops import All, Measure, QubitOperator, TimeEvolution, X, H


import projectq.setups.decompositions as rules
from projectq.cengines import (AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               DecompositionRuleSet)
# ====================
# from hiq.projectq.cengines import GreedyScheduler, HiQMainEngine
# from hiq.projectq.backends import SimulatorMPI

# ====================

from projectq.backends import Simulator
from projectq import MainEngine

from mpi4py import MPI

import datetime



def H_Ising(N, J, alpha, beta):
    '''
    implement the target Hamiltonian given in the probelm
    parameters
        N: (int) number of allocated qubits
        J, alpha, beta: (float) coupling coefficients
    return
        H:(QubitOperator) target hamiltonian
    '''
    H = 0*QubitOperator("")
    for i in range(N):
        H += alpha*QubitOperator("X"+str(i))
        H += beta*QubitOperator("Z"+str(i))
    for i in range(N-1):
        H += J*QubitOperator("X{0} X{1}".format(i, i+1))
        H += J*QubitOperator("Z{0} Z{1}".format(i, i+1))
    return (-1)*H


def H_ansatz(N, B):
    '''
    implement the mixer Hamiltonian given in the probelm
    parameters
        N: (int) number of allocated qubits
        B: (float) coupling coefficients
    return
        H:(QubitOperator) target hamiltonian
    '''
    H = 0*QubitOperator("")
    for i in range(N):
        H += B*QubitOperator("X"+str(i))
    return (-1)*H



# contruct circuits to implement discreted annealing
if __name__ == "__main__":
    # define parameters
    N=13 # qubit number 
    M=300 # iteration 
    J=2 
    alpha=0.5
    beta=1.0
    B=0.5
    T = 6/J
    dt = T/M

    # initialize engine

    backend = Simulator(gate_fusion=True)

    cache_depth = 10
    rule_set = DecompositionRuleSet(modules=[rules])
    engines = [TagRemover()
                , LocalOptimizer(cache_depth)
                , AutoReplacer(rule_set)
                , TagRemover()
                , LocalOptimizer(cache_depth)
                ]

    engine = MainEngine(backend, engines)

    # backend = SimulatorMPI(gate_fusion=True, num_local_qubits=N)

    # cache_depth = M*3
    # rule_set = DecompositionRuleSet(modules=[rules])
    # engines = [TagRemover(), LocalOptimizer(cache_depth), AutoReplacer(rule_set), TagRemover(), LocalOptimizer(cache_depth), GreedyScheduler()
    #         ]

    # engine = HiQMainEngine(backend, engines)
    start=datetime.datetime.now()
    print(start)


    print("initialize engine..")

    # define energy hamiltonian
    H_t=H_Ising(N,J,alpha,beta)
    H_0=H_ansatz(N,B)

    # start annealing 
    print("start annealing..")

    state = engine.allocate_qureg(N)
    All(H) | state  # |+]



    j=0
    for i in range(N):
        TimeEvolution(-B*dt/2, QubitOperator("X"+str(i))) | state

    # j= 1:M-1
    for j in range(1, M):
        for i in range(N):
            TimeEvolution(-beta*(j*dt**2/T), QubitOperator("Z"+str(i))) | state
            TimeEvolution(-(alpha*(j*dt**2/T)+B*dt*(1-j*dt/T)),QubitOperator("X"+str(i))) | state
        for i in range(N-1):
            TimeEvolution(-J*(j*dt**2/T), QubitOperator("Z{0} Z{1}".format(i, i+1))) | state
            TimeEvolution(-J*(j*dt**2/T), QubitOperator("X{0} X{1}".format(i, i+1))) | state
    # j=M
    for i in range(N):
        TimeEvolution(-beta*(dt/2), QubitOperator("Z"+str(i))) | state
        TimeEvolution(-alpha*(dt/2), QubitOperator("X"+str(i))) | state
    for i in range(N-1):
        TimeEvolution(-J*(dt/2), QubitOperator("Z{0} Z{1}".format(i, i+1))) | state
        TimeEvolution(-J*(dt/2), QubitOperator("X{0} X{1}".format(i, i+1))) | state

    engine.flush()

    end=datetime.datetime.now()
    print(end)
    print("time used:",(end-start))
    energy=engine.backend.get_expectation_value(H_t,state)
    All(Measure)|state
    print("final energy:",energy)
