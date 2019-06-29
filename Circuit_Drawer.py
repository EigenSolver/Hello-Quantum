from Bell_States_Generator import bell_state
from Quantum_Teleportation import quantum_teleportation
from projectq import MainEngine
from projectq.backends import CircuitDrawer
from projectq.ops import H,Rz
import os

drawer=CircuitDrawer()
eng=MainEngine(drawer)
qubit_psi=eng.allocate_qubit()
H|qubit_psi
Rz(1.32)|qubit_psi
quantum_teleportation(eng,qubit_psi)
eng.flush()

with open("./circuit_plot/quantum_teleportation.tex","w") as f:
    f.write(drawer.get_latex())