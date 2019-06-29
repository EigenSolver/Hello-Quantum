from projectq import MainEngine
from projectq.ops import H,Measure


def quantum_random_number_generator(eng): 
    qubit=eng.allocate_qubit()

    # apply Hadamard gate and generate superposition state
    H|qubit
    # Measure the qubit
    Measure|qubit

    return int(qubit)


if __name__=="__main__":
     eng=MainEngine()
     measurements=[]

     N=10000

     for i in range(N):    
          measurements.append(quantum_random_number_generator(eng))
     # end engine
     eng.flush()

     ones=sum(measurements)
     zeros=N-ones

     print("By allocating {} qubits and applying Hadamard gates on them,\
          we do measurements on the superpositions and get the ratio between 1 and 0: ".format(N))
     print(ones,":",zeros)
