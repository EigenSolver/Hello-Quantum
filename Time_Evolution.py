from projectq.ops import All, QubitOperator, H, Measure, X, Z, TimeEvolution
from projectq import MainEngine
import numpy as np

eng=MainEngine()
state=eng.allocate_qureg(2)

All(H)|state

# TimeEvolution(np.pi/2,QubitOperator("X0")+QubitOperator("Z1"))|state
# it can be prove that the two approach above is equal

TimeEvolution(np.pi/2,QubitOperator("X0"))|state
TimeEvolution(np.pi/2,QubitOperator("Z1"))|state
eng.flush()
print(eng.backend.cheat())
All(Measure)|state


# print([int(x) for x in state])
