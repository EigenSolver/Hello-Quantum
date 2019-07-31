# -*- coding: utf-8 -*- 
# File name: problem1.py

# Huawei Quantum Programming Competition
# Author: QuantumSupremacy Group

import numpy as np
import datetime
import projectq as pq
from projectq import MainEngine
from projectq.ops import All, Measure, QubitOperator, TimeEvolution, X, H

import projectq.setups.decompositions as rules
from projectq.cengines import (AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               DecompositionRuleSet)

from projectq.backends import Simulator
from hiq.projectq.backends import SimulatorMPI
from hiq.projectq.cengines import GreedyScheduler, HiQMainEngine

from mpi4py import MPI


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
    H = (-1)*H
    return H

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
        H += QubitOperator("X"+str(i))
    return (-B)*H

def annealing(engine,M):
    '''
    annealing with Trotter coefficient M
    Returns:
        energy(float): energy expectation of end state
    '''
    # define parameters
    N=13 # qubit number 
    #M: iteration 
    J=2 
    alpha=0.5
    beta=1.0
    B=0.5
    T = 6/J
    dt = T/M
    
    H_t=H_Ising(N,J,alpha,beta)
    H_0=H_ansatz(N,B)

    state = engine.allocate_qureg(N)
    All(H) | state  # |+]
    
    TimeEvolution(dt/2,H_0)|state
    for j in range(1,M-1):
        TimeEvolution(j*dt**2/T,H_t)|state
        TimeEvolution(dt*(1-j*dt/T),H_0)|state
    TimeEvolution(dt/2, H_t)|state


    engine.flush()

    energy=engine.backend.get_expectation_value(H_t,state)
    All(Measure)|state
    
    return energy
    
def adiabatic_simulation(engine):
    """
    The function you need to modify. Returns:

    real_energy(float):

    The final ideally continously evolved energy.

    simulated_energy(float):

    The final energy simulated by your model.
    """
    
    simulated_energy=annealing(engine,20)
    real_energy = annealing(engine,100)
    return simulated_energy, real_energy


if __name__ == "__main__":
    # use projectq simulator
    #eng = MainEngine()

    # use hiq simulator
    # backend = SimulatorMPI(gate_fusion=True, num_local_qubits=13)

    cache_depth = 10
    rule_set = DecompositionRuleSet(modules=[rules])
    engines = [TagRemover(), LocalOptimizer(cache_depth), AutoReplacer(rule_set), TagRemover(), LocalOptimizer(cache_depth)]

    # make the compiler and run the circuit on the simulator backend
    backend = Simulator(gate_fusion=True)
    eng = HiQMainEngine(backend, engines)

    simulated_energy, real_energy = adiabatic_simulation(eng)

    simulated_error = simulated_energy - real_energy

    print(simulated_error/real_energy)
