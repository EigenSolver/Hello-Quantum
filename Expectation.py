import projectq as pq
from projectq import MainEngine
from projectq.ops import All
from projectq.ops import H, X, Z, Measure, QubitOperator
import numpy as np

eng=MainEngine()
qb_array=eng.allocate_qureg(4)

All(H)|qb_array

eng.flush()

operator=QubitOperator("X0 X1")
expectation=eng.backend.get_expectation_value(operator,qb_array)
All(Measure)|qb_array
print("X_0*X_1",expectation)

# ===============================

qb_array=eng.allocate_qureg(3)

All(H)|qb_array
All(Measure)|qb_array
eng.flush()

print([int(qb) for qb in qb_array])

operator=QubitOperator("Z0")+QubitOperator("Z1")+QubitOperator("Z2")
expectation=eng.backend.get_expectation_value(operator,qb_array)

All(Measure)|qb_array
print("Z_0+Z_1+Z_2",expectation)
