from projectq.ops import All, NOT, CNOT, H, Measure, Rz, X, Z
from projectq import MainEngine


def bell_state(eng,inputs=(0,0)):
    bit1=eng.allocate_qubit()
    bit2=eng.allocate_qubit()
    

    if inputs[0]:
        NOT|bit1
    if inputs[1]:
        NOT|bit2

    H|bit1
    CNOT|(bit1,bit2)

    return bit1,bit2

if __name__=="__main__":
    eng=MainEngine()

    input_states=[(0,0),(0,1),(1,0),(1,1)]
    # generate different bell states with corresponding input 
    for input_state in input_states:
        print("With given input states {}, \
                we measure the result as following:".format(input_state))
        
        # simple test to show the entanglement of the Bell state
        for i in range(10):
            (b1,b2)=bell_state(eng,input_state)

            Measure|b1
            Measure|b2
            print(int(b1),int(b2))

    eng.flush()
