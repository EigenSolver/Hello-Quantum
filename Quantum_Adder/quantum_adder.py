from projectq.ops import H, X, Z, Rx, Ry, Rz, Measure, QubitOperator, C, CNOT, All
from projectq.meta import Loop, Compute, Uncompute, Control
import projectq
import numpy as np

#=================================

import projectq.setups.decompositions as rules

from projectq.backends import Simulator

from projectq.cengines import (AutoReplacer,
                            LocalOptimizer,
                            TagRemover,
                            DecompositionRuleSet)

from projectq import MainEngine

backend = Simulator(gate_fusion=True)

cache_depth = 10
rule_set = DecompositionRuleSet(modules=[rules])
engines = [TagRemover(), LocalOptimizer(cache_depth), AutoReplacer(rule_set), TagRemover(), LocalOptimizer(cache_depth)
        ]

engine = MainEngine(backend, engines)

# ================================

def quantum_fourier_transform(state):
    '''
    parameters
        state: (qureg) quantum register representing a quantum state
    '''
    assert type(state)==projectq.types.Qureg

    n=len(state)

    for i in range(n):
        H|state[i]
        for j in range(i+1,n):
            C(Rz(np.pi/(2**j)))|(state[j],state[i])

def ancillary_add(trans_state, ori_state):
    '''
    parameters
        trans_state: (qureg) state a after QFT
        ori_state: (qureg) original state b, allocated in the same engine
    '''
    assert len(trans_state)==len(ori_state)
    n=len(trans_state)

    for i in range(n):
        for j in range(i,n):
            C(Rz(np.pi/(2**j)))|(ori_state[j],trans_state[i])

def Quantum_Adder(state_a,state_b,eng=engine):
    '''
    parameters
        state_a:(qureg)
        state_b:(qureg)
    '''
    with Compute(eng):
        quantum_fourier_transform(state_a)
    ancillary_add(state_a,state_b)
    Uncompute(eng)
    # return state_a

def prep_state(bitstring,eng=engine):
    n=len(bitstring)

    state=eng.allocate_qureg(n)
    for i in range(n):
        if bitstring[i]:
            X|state[i]
    return state

# ========== test ================
if __name__=="__main__":

    a=[1,1,0,0,1]
    print("a=",a)
    state_a=prep_state(a)


    b=[0,0,1,0,1]
    print("b=",b)
    state_b=prep_state(b)

    Quantum_Adder(state_a,state_b)

    All(Measure)|state_a
    engine.flush()
    print("quantum adding...")
    print("a+b=",[int(qb) for qb in state_a])

    

    
