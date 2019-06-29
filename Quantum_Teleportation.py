from projectq.ops import All, CNOT, H, Measure, Rz, X, Z
from projectq import MainEngine
from projectq.meta import Dagger, Control

from Bell_States_Generator import bell_state
from random import random

def quantum_teleportation(eng,qubit_psi):
    # generate entangled state betweeen Alice and Bob
    qubit_A, qubit_B=bell_state(eng)

    # suppose we have a unknown gerenate function to produce the psi state 
    # which is to be teleported 
    
    CNOT|(qubit_psi,qubit_A) # notice that the first input bit is the control bit
    # Hadamard
    H|qubit_psi
    Measure|qubit_psi
    Measure|qubit_A

    with Control(eng,qubit_A):
        X|qubit_B
    with Control(eng,qubit_psi):
        Z|qubit_B
    
    return qubit_B

def random_qubit(eng):
    qb=eng.allocate_qubit()
    H|qb
    Rz(random())|qb
    return qb


if __name__=="__main__":
    def f(eng,qubit):
        H|qubit
        Rz(0.618)|qubit

    eng=MainEngine()
    psi=eng.allocate_qubit()
    f(eng,psi)

    teleported_state=quantum_teleportation(eng,psi)
    # print(teleported_state==psi)

    with Dagger(eng):
        f(eng,teleported_state)
    del teleported_state

    Measure|psi
    eng.flush()
    print("succeed!")



