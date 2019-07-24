import projectq as pq
from projectq import MainEngine
from projectq.ops import All
from projectq.ops import H, X, Z, Measure, QubitOperator
import numpy as np

n = 200
count=[]
for i in range(n):
    eng = MainEngine()

    psi = eng.allocate_qubit()
    H | psi
    Measure | psi
    eng.flush()
    count.append(int(psi))

print("expectations: {}".format(np.mean(count)))


configs=[]
for i in range(n):
    eng = MainEngine()
    wave_function = eng.allocate_qureg(2)
    All(H) | wave_function

    # #???
    All(Measure)|wave_function
    eng.flush()
    configs.append([int(qb) for qb in wave_function])

print("10: ",configs.count([1,0]),"/",n)
print("11: ",configs.count([1,1]),"/",n)
print("00: ",configs.count([0,0]),"/",n)
print("01: ",configs.count([0,1]),"/",n)
# print("expectation:", energy)

