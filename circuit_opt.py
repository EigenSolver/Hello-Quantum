import time

import numpy as np

from hiq.backends import CircuitDrawer, CommandPrinter, Simulator

from hiq.cengines import (MainEngine,

                               AutoReplacer,

                               LocalOptimizer,

                               TagRemover,

                               DecompositionRuleSet)

 

from huawei.hiq.cengines import GreedyScheduler

import hiq.setups.decompositions

 

import importlib

import hiq.setups.default

from hiq.meta import Dagger, Control

from hiq.ops import H, Z, X, Measure, All, C, T, SqrtX, SqrtY

from hiq.libs.math.entropy import entropy

from hiq.ops import All, H

 

import ctypes

ctypes.CDLL("libmpi.so", mode=ctypes.RTLD_GLOBAL)

 

from huawei.hiq.cengines import GreedyScheduler

 

nqubits = 25

d = 30

def run(eng):

    # Alice creates a nice state to send

    x = eng.allocate_qureg(nqubits)
    # Start with a cycle of Hadamard gates

    All(H) | x
    
    # Cycle #0

    C(Z) | (x[0], x[5])

    C(Z) | (x[11], x[16])

    C(Z) | (x[2], x[7])

    C(Z) | (x[13], x[18])

    C(Z) | (x[4], x[9])
    
    
    # Cycle #1

    C(Z) | (x[10], x[15])

    C(Z) | (x[1], x[6])

    C(Z) | (x[12], x[17])

    C(Z) | (x[3], x[8])

    C(Z) | (x[14], x[19])

    SqrtY | x[0]

    T | x[5]

    SqrtX | x[11]

    T | x[16]

    SqrtX | x[2]

    SqrtY | x[7]

    SqrtX | x[13]

    SqrtX | x[18]

    T | x[4]

    SqrtX | x[9]
    
    # Cycle #2

    C(Z) | (x[5], x[10])

    C(Z) | (x[16], x[21])

    C(Z) | (x[7], x[12])

    C(Z) | (x[18], x[23])

    C(Z) | (x[9], x[14])

    T | x[15]

    SqrtX | x[1]

    T | x[6]

    T | x[17]

    SqrtY | x[3]

    T | x[8]

    SqrtY | x[19]
    
    # Cycle #3

    C(Z) | (x[15], x[20])

    C(Z) | (x[6], x[11])

    C(Z) | (x[17], x[22])

    C(Z) | (x[8], x[13])

    C(Z) | (x[19], x[24])

    SqrtY | x[5]

    SqrtX | x[10]

    SqrtX | x[16]

    SqrtY | x[21]

    SqrtY | x[7]

    SqrtY | x[12]

    SqrtX | x[18]

    T | x[23]

    SqrtX | x[9]

    SqrtY | x[14]
    
    # Cycle #4

    C(Z) | (x[0], x[1])

    C(Z) | (x[7], x[8])

    C(Z) | (x[10], x[11])

    C(Z) | (x[17], x[18])

    C(Z) | (x[20], x[21])

    SqrtX | x[15]

    SqrtY | x[6]

    T | x[22]

    SqrtY | x[13]

    SqrtY | x[19]

    SqrtX | x[24]
    
    # Cycle #5

    C(Z) | (x[2], x[3])

    C(Z) | (x[5], x[6])

    C(Z) | (x[12], x[13])

    C(Z) | (x[15], x[16])

    C(Z) | (x[22], x[23])

    SqrtX | x[0]

    SqrtY | x[1]

    SqrtX | x[7]

    SqrtY | x[8]

    SqrtY | x[10]

    SqrtX | x[11]

    SqrtX | x[17]

    SqrtX | x[18]

    T | x[20]

    SqrtY | x[21]
    
    # Cycle #6

    C(Z) | (x[1], x[2])

    C(Z) | (x[8], x[9])

    C(Z) | (x[11], x[12])

    C(Z) | (x[18], x[19])

    C(Z) | (x[21], x[22])

    SqrtX | x[3]

    SqrtX | x[5]

    SqrtY | x[6]

    SqrtY | x[13]

    SqrtY | x[15]

    SqrtY | x[16]

    SqrtY | x[23]
    
    # Cycle #7

    C(Z) | (x[3], x[4])

    C(Z) | (x[6], x[7])

    C(Z) | (x[13], x[14])

    C(Z) | (x[16], x[17])

    C(Z) | (x[23], x[24])

    SqrtX | x[1]

    SqrtY | x[2]

    SqrtX | x[8]

    SqrtY | x[9]

    SqrtY | x[11]

    SqrtY | x[12]

    SqrtY | x[18]

    SqrtX | x[19]

    SqrtX | x[21]

    SqrtY | x[22]
    
    # Cycle #8

    C(Z) | (x[0], x[5])

    C(Z) | (x[11], x[16])

    C(Z) | (x[2], x[7])

    C(Z) | (x[13], x[18])

    C(Z) | (x[4], x[9])

    SqrtY | x[3]

    SqrtX | x[6]

    SqrtX | x[14]

    SqrtX | x[17]

    SqrtY | x[23]

    SqrtY | x[24]

    # Cycle #9

    C(Z) | (x[10], x[15])

    C(Z) | (x[1], x[6])

    C(Z) | (x[12], x[17])

    C(Z) | (x[3], x[8])

    C(Z) | (x[14], x[19])

    SqrtY | x[0]

    SqrtY | x[5]

    SqrtY | x[11]

    SqrtY | x[16]

    SqrtX | x[2]

    SqrtX | x[7]

    SqrtY | x[13]

    SqrtX | x[18]

    SqrtX | x[4]

    SqrtY | x[9]
    
    
    # Cycle #10

    C(Z) | (x[5], x[10])

    C(Z) | (x[16], x[21])

    C(Z) | (x[7], x[12])

    C(Z) | (x[18], x[23])

    C(Z) | (x[9], x[14])

    SqrtX | x[15]

    SqrtY | x[1]

    SqrtY | x[6]

    SqrtX | x[17]

    SqrtX | x[3]

    SqrtY | x[8]

    SqrtY | x[19]

    # Cycle #11

    C(Z) | (x[15], x[20])

    C(Z) | (x[6], x[11])

    C(Z) | (x[17], x[22])

    C(Z) | (x[8], x[13])

    C(Z) | (x[19], x[24])

    SqrtX | x[5]

    SqrtX | x[10]

    SqrtX | x[16]

    SqrtX | x[21]

    SqrtX | x[7]

    SqrtX | x[12]

    SqrtY | x[18]

    SqrtY | x[23]

    SqrtX | x[9]

    SqrtY | x[14]
    
    # Cycle #12

    C(Z) | (x[0], x[1])

    C(Z) | (x[7], x[8])

    C(Z) | (x[10], x[11])

    C(Z) | (x[17], x[18])

    C(Z) | (x[20], x[21])

    SqrtY | x[15]

    SqrtY | x[6]

    SqrtX | x[22]

    SqrtY | x[13]

    SqrtX | x[19]

    SqrtX | x[24]
    
    # Cycle #13

    C(Z) | (x[2], x[3])

    C(Z) | (x[5], x[6])

    C(Z) | (x[12], x[13])

    C(Z) | (x[15], x[16])

    C(Z) | (x[22], x[23])

    SqrtX | x[0]

    SqrtY | x[1]

    SqrtY | x[7]

    SqrtY | x[8]

    SqrtY | x[10]

    SqrtX | x[11]

    SqrtY | x[17]

    SqrtX | x[18]

    SqrtY | x[20]

    SqrtX | x[21]
    
    # Cycle #14

    C(Z) | (x[1], x[2])

    C(Z) | (x[8], x[9])

    C(Z) | (x[11], x[12])

    C(Z) | (x[18], x[19])

    C(Z) | (x[21], x[22])

    SqrtY | x[3]

    SqrtX | x[5]

    SqrtX | x[6]

    SqrtY | x[13]

    SqrtY | x[15]

    SqrtY | x[16]

    SqrtX | x[23]

    # Cycle #15

    C(Z) | (x[3], x[4])

    C(Z) | (x[6], x[7])

    C(Z) | (x[13], x[14])

    C(Z) | (x[16], x[17])

    C(Z) | (x[23], x[24])

    SqrtX | x[1]

    SqrtY | x[2]

    SqrtY | x[8]

    SqrtX | x[9]

    SqrtY | x[11]

    SqrtY | x[12]

    SqrtX | x[18]

    SqrtY | x[19]

    SqrtX | x[21]

    SqrtX | x[22]

 

    # Cycle #16

    C(Z) | (x[0], x[5])

    C(Z) | (x[11], x[16])

    C(Z) | (x[2], x[7])

    C(Z) | (x[13], x[18])

    C(Z) | (x[4], x[9])

    SqrtY | x[3]

    SqrtX | x[6]

    SqrtX | x[14]

    SqrtX | x[17]

    SqrtX | x[23]

    SqrtY | x[24]

 

    # Cycle #17

    C(Z) | (x[10], x[15])

    C(Z) | (x[1], x[6])

    C(Z) | (x[12], x[17])

    C(Z) | (x[3], x[8])

    C(Z) | (x[14], x[19])

    SqrtY | x[0]

    SqrtY | x[5]

    SqrtY | x[11]

    SqrtY | x[16]

    SqrtX | x[2]

    SqrtY | x[7]

    SqrtX | x[13]

    SqrtY | x[18]

    SqrtX | x[4]

    SqrtY | x[9]

 

    # Cycle #18

    C(Z) | (x[5], x[10])

    C(Z) | (x[16], x[21])

    C(Z) | (x[7], x[12])

    C(Z) | (x[18], x[23])

    C(Z) | (x[9], x[14])

    SqrtY | x[15]

    SqrtY | x[1]

    SqrtY | x[6]

    SqrtY | x[17]

    SqrtY | x[3]

    SqrtX | x[8]

    SqrtY | x[19]

 

    # Cycle #19

    C(Z) | (x[15], x[20])

    C(Z) | (x[6], x[11])

    C(Z) | (x[17], x[22])

    C(Z) | (x[8], x[13])

    C(Z) | (x[19], x[24])

    SqrtX | x[5]

    SqrtY | x[10]

    SqrtY | x[16]

    SqrtY | x[21]

    SqrtX | x[7]

    SqrtX | x[12]

    SqrtY | x[18]

    SqrtX | x[23]

    SqrtX | x[9]

    SqrtX | x[14]

 

    # Cycle #20

    C(Z) | (x[0], x[1])

    C(Z) | (x[7], x[8])

    C(Z) | (x[10], x[11])

    C(Z) | (x[17], x[18])

    C(Z) | (x[20], x[21])

    SqrtY | x[15]

    SqrtY | x[6]

    SqrtX | x[22]

    SqrtY | x[13]

    SqrtX | x[19]

    SqrtX | x[24]

 

    # Cycle #21

    C(Z) | (x[2], x[3])

    C(Z) | (x[5], x[6])

    C(Z) | (x[12], x[13])

    C(Z) | (x[15], x[16])

    C(Z) | (x[22], x[23])

    SqrtX | x[0]

    SqrtX | x[1]

    SqrtY | x[7]

    SqrtX | x[8]

    SqrtY | x[10]

    SqrtY | x[11]

    SqrtX | x[17]

    SqrtX | x[18]

    SqrtY | x[20]

    SqrtY | x[21]

 

    # Cycle #22

    C(Z) | (x[1], x[2])

    C(Z) | (x[8], x[9])

    C(Z) | (x[11], x[12])

    C(Z) | (x[18], x[19])

    C(Z) | (x[21], x[22])

    SqrtX | x[3]

    SqrtY | x[5]

    SqrtY | x[6]

    SqrtX | x[13]

    SqrtX | x[15]

    SqrtX | x[16]

    SqrtX | x[23]

 

    # Cycle #23

    C(Z) | (x[3], x[4])

    C(Z) | (x[6], x[7])

    C(Z) | (x[13], x[14])

    C(Z) | (x[16], x[17])

    C(Z) | (x[23], x[24])

    SqrtY | x[1]

    SqrtX | x[2]

    SqrtX | x[8]

    SqrtX | x[9]

    SqrtY | x[11]

    SqrtX | x[12]

    SqrtX | x[18]

    SqrtY | x[19]

    SqrtY | x[21]

    SqrtY | x[22]

 

    # Cycle #24

    C(Z) | (x[0], x[5])

    C(Z) | (x[11], x[16])

    C(Z) | (x[2], x[7])

    C(Z) | (x[13], x[18])

    C(Z) | (x[4], x[9])

    SqrtY | x[3]

    SqrtX | x[6]

    SqrtY | x[14]

    SqrtX | x[17]

    SqrtX | x[23]

    SqrtX | x[24]

 

    # Cycle #25

    C(Z) | (x[10], x[15])

    C(Z) | (x[1], x[6])

    C(Z) | (x[12], x[17])

    C(Z) | (x[3], x[8])

    C(Z) | (x[14], x[19])

    SqrtY | x[0]

    SqrtY | x[5]

    SqrtY | x[11]

    SqrtX | x[16]

    SqrtX | x[2]

    SqrtX | x[7]

    SqrtX | x[13]

    SqrtX | x[18]

    SqrtX | x[4]

    SqrtX | x[9]

 

    # Cycle #26

    C(Z) | (x[5], x[10])

    C(Z) | (x[16], x[21])

    C(Z) | (x[7], x[12])

    C(Z) | (x[18], x[23])

    C(Z) | (x[9], x[14])

    SqrtX | x[15]

    SqrtX | x[1]

    SqrtX | x[6]

    SqrtX | x[17]

    SqrtY | x[3]

    SqrtX | x[8]

    SqrtX | x[19]

 

    # Cycle #27

    C(Z) | (x[15], x[20])

    C(Z) | (x[6], x[11])

    C(Z) | (x[17], x[22])

    C(Z) | (x[8], x[13])

    C(Z) | (x[19], x[24])

    SqrtY | x[5]

    SqrtY | x[10]

    SqrtX | x[16]

    SqrtY | x[21]

    SqrtY | x[7]

    SqrtY | x[12]

    SqrtX | x[18]

    SqrtX | x[23]

    SqrtY | x[9]

    SqrtX | x[14]

 

    # Cycle #28

    C(Z) | (x[0], x[1])

    C(Z) | (x[7], x[8])

    C(Z) | (x[10], x[11])

    C(Z) | (x[17], x[18])

    C(Z) | (x[20], x[21])

    SqrtY | x[15]

    SqrtY | x[6]

    SqrtX | x[22]

    SqrtY | x[13]

    SqrtX | x[19]

    SqrtY | x[24]
    # Cycle #29

    C(Z) | (x[2], x[3])

    C(Z) | (x[5], x[6])

    C(Z) | (x[12], x[13])

    C(Z) | (x[15], x[16])

    C(Z) | (x[22], x[23])

    SqrtY | x[0]

    SqrtX | x[1]

    SqrtY | x[7]

    SqrtY | x[8]

    SqrtX | x[10]

    SqrtX | x[11]

    SqrtY | x[17]

    SqrtY | x[18]

    SqrtX | x[20]
    
    SqrtY | x[21]

drawer=CircuitDrawer()
eng=MainEngine(drawer)

run(eng)
with open("./circuit_plot/circuit_opt.tex","w") as f:
    f.write(drawer.get_latex())

