from projectq.ops import H, X, Z, Rx, Ry, Rz, Measure, QubitOperator, CNOT, All
from projectq.meta import Loop, Compute, Uncompute, Control

import numpy as np

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

# ============================ header


def bitstring(n):
    if n > 1:
        return ['0'+bits for bits in bitstring(n-1)]+['1'+bits for bits in bitstring(n-1)]
    elif n == 1:
        return ['0', '1']
    else:
        print('error')

def cost(param, target_state, p, eng):
    '''
    VQE core part
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

        CNOT | (qureg[n-1], qureg[0])
        for i in range(n-1):
            CNOT | (qureg[i], qureg[i+1])

    eng.flush()
    for i in range(len(target_state)):
        count += target_state[i] * eng.backend.get_amplitude(all_string[i], qureg)
    All(Measure) | qureg

    return 1-np.abs(count)
# ===============================
def check_solution(param,p,eng):
    n=int(len(param)/3/p)
    param=param.reshape((n,3*p))
    amplitudes=[]

    qureg = eng.allocate_qureg(n)
    All(H) | qureg

    for m in range(p):
        for i in range(n):
            Rx(param[i, 3*m+0]) | qureg[i]
            Ry(param[i, 3*m+1]) | qureg[i]
            Rz(param[i, 3*m+2]) | qureg[i]

        CNOT | (qureg[n-1], qureg[0])
        for i in range(n-1):
            CNOT | (qureg[i], qureg[i+1])

    eng.flush()
    amplitudes=[eng.backend.get_amplitude(conf, qureg) for conf in bitstring(n)]

    All(Measure)|qureg
    return amplitudes

# ====================================

from scipy.optimize import minimize

def VQC(target_state,p, eng, method="COBYLA"):
    p = 3
    n = int(np.log2(len(target_state)))
    param = np.ones((n, 3*p)).flatten()
    result=minimize(cost,x0=param, args=(target_state,p,engine),method=method)#"TNC")#
    return result
        

