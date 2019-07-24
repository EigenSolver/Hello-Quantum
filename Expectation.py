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

print(expectation)