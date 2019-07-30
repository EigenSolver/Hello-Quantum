import projectq as pq
from projectq import MainEngine
from projectq.ops import All
from projectq.ops import H, X, Z, Measure, QubitOperator
import numpy as np

eng=MainEngine()
qb_array=eng.allocate_qubit()

# X|qb_array

eng.flush()
exp=eng.backend.get_expectation_value(Z,qb_array)
Measure|qb_array
print([int(q) for q in qb_array])
print(exp)