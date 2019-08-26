# file name: problem1.py
from projectq import MainEngine
from projectq.ops import X, H, CX, All, Measure, CX
import numpy as np

def subproblem1():
    eng = MainEngine()
    # Allocate qubits
    qubits = eng.allocate_qureg(2)
    # You should add quantum operations here to change the state of qubits into |01>
    # X|qubits[1]
    X|qubits[0]
    # Submit all your operations
    eng.flush()
    # Get the amplitudes list
    amplitudes = np.array(eng.backend.cheat()[1])
    print("p1:{}".format(amplitudes))
    # Here we return the amplitudes without measures and you do not need to care about the warnings
    All(Measure) | qubits
    return amplitudes

def subproblem2():
    eng = MainEngine()
    qubits = eng.allocate_qureg(2)
    # You should add quantum operations here to change the state of qubits into |00>+|11>
    H|qubits[0]
    CX|(qubits[0],qubits[1])

    eng.flush()
    amplitudes = np.array(eng.backend.cheat()[1])
    print("p2:{}".format(amplitudes))
    All(Measure) | qubits
    return amplitudes

def subproblem3():
    eng = MainEngine()
    qubits = eng.allocate_qureg(2)
    # You should add quantum operations here to change the state of qubits into |10>+|01>
    X|qubits[1]

    H|qubits[0]
    CX|(qubits[0],qubits[1])

    eng.flush()
    amplitudes = np.array(eng.backend.cheat()[1])
    print("p3:{}".format(amplitudes))
    All(Measure) | qubits
    return amplitudes

def subproblem4():
    eng = MainEngine()
    qubits = eng.allocate_qureg(3)
    # You should add quantum operations here to change the state of qubits into |101>+|010>
    H|qubits[1]

    CX|(qubits[1],qubits[0])
    CX|(qubits[1],qubits[2])
    
    eng.flush()
    amplitudes = np.array(eng.backend.cheat()[1])
    All(Measure) | qubits
    print("p4:{}".format(amplitudes))
    return amplitudes

if __name__ == '__main__':
    # You can use the main function to evaluate your own score.
    # Initialize your score
    score = 0
    amplitude1 = subproblem1()
    target_amplitude1 = np.array([0, 1, 0, 0])
    if sum(np.abs(amplitude1 - target_amplitude1)) <= 0.001:
        score += 5
    amplitude2 = subproblem2()
    target_amplitude2 = np.array([np.sqrt(2)/2, 0, 0, np.sqrt(2)/2])
    if sum(np.abs(amplitude2 - target_amplitude2)) <= 0.001:
        score += 5
    amplitude3 = subproblem3()
    target_amplitude3 = np.array([0, np.sqrt(2)/2, np.sqrt(2)/2, 0])
    if sum(np.abs(amplitude3 - target_amplitude3)) <= 0.001:
        score += 5
    amplitude4 = subproblem4()
    target_amplitude4 = np.array([0, 0, np.sqrt(2)/2, 0, 0, np.sqrt(2)/2, 0, 0])
    if sum(np.abs(amplitude4 - target_amplitude4)) <= 0.001:
        score += 5
    print ('Your score of problem1 is:{}'.format(score))
