from Variational_Quantum_Circuit import *
import numpy as np
import pandas as pd 

N=30
n=4
p=3
method="COBYLA"

def normalize(array):
    return array/np.sqrt(np.sum(array**2))

def gen_random_state(n):
    return normalize(np.array([np.random.rand() for i in range(n)]))


# save the following data
# target state
# generated state
# distance
# optimal parameter
# iteration


# engine is defined in VQC
def report(i,N,report_rate=2):
    if (i+1)%report_rate==0:
        print("progress {0}/{1}".format(i+1,N))

data=[]
print("starting benchmarking...")
for i in range(N):
    target=gen_random_state(n)
    result=VQC(target,p,engine,method=method)
    cost=result.fun
    optimal_param=result.x
    n_iter=result.nfev
    gen_state=check_solution(optimal_param,p,engine)
    data.append([target,gen_state,optimal_param,cost,n_iter])
    report(i,N,5)

print("writing data...")

wrapper=pd.DataFrame(data,columns=["target","result","opt_param","cost","n_iter"])

wrapper.to_csv("benchmark_n={0}_p={1}_method={2}.csv".format(n,p,method))