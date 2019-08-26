# file name: problem2.py
def example1():
    """ This is an example for matrix [[0, 1], [1, 0]]
    Note:
        For the multi-qubit or multi-operations quantum circuit, use ';' to split them.
    """
    quantum_circuit = 'X | qubits[0]'
    return quantum_circuit

def example2():
    """ This is an example for matrix [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
    Note:
        For the multi-qubit or multi-operations quantum circuit, use ';' to split them.
    """
    quantum_circuit = 'CX | (qubits[0], qubits[1])'
    return quantum_circuit

def example3():
    """ This is an example for matrix [[0, 0, 0, 1], [0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0]]
    Note:
        For the multi-qubit or multi-operations quantum circuit, use ';' to split them.
    """
    quantum_circuit = 'X | qubits[0]; X | qubits[1]; CX | (qubits[0], qubits[1])'
    return quantum_circuit

def subproblem1():
    """ Now you need to find a circuit equals to matrix [[1, 1], [1, -1]] / sqrt(2)
    """
    # You can change your circuit here
    quantum_circuit = 'H | qubits[0]'
    return quantum_circuit

def subproblem2():
    """ Now you need to find a circuit equals to matrix [[1, 1], [-1, 1]] / sqrt(2)
    """
    # You can change your circuit here
    quantum_circuit = 'X | qubits[0]; H | qubits[0]'
    return quantum_circuit

def subproblem3():
    """ Now you need to find a circuit equals to matrix [[0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0], [1, 0, 0, 0]]
    """
    # You can change your circuit here
    quantum_circuit = 'X | qubits[0]; X | qubits[1])'
    return quantum_circuit

def subproblem4():
    """ Now you need to find a circuit equals to matrix [[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]]
    """
    # You can change your circuit here
    quantum_circuit = 'CX | (qubits[0], qubits[1]); CX | (qubits[1], qubits[0]); CX | (qubits[0], qubits[1])'
    return quantum_circuit

def subproblem5():
    """ Now you need to find a circuit equals to matrix [[1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0], [0, 1, 0, 0]]
    """
    # You can change your circuit here
    quantum_circuit = ''
    return quantum_circuit
