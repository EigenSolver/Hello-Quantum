# -*- coding: utf-8 -*-
# File name: problem3.py

# Huawei Quantum Programming Competition
# Author: QuantumSupremacy Group
# Github: https://github.com/Neuromancer43/Hello-Quantum/tree/master/HiQ_competition
from __future__ import division
import random


import numpy as np
from scipy.optimize import minimize

import projectq as pq
from projectq import MainEngine
from projectq.ops import X, Y, Z, H, S, T, CX, CZ, Rx, Ry, Rz, Measure, All


import projectq.setups.decompositions as rules
from projectq.cengines import (AutoReplacer,
                               LocalOptimizer,
                               TagRemover,
                               DecompositionRuleSet)
# ====================
from hiq.projectq.cengines import GreedyScheduler, HiQMainEngine
from hiq.projectq.backends import SimulatorMPI

from mpi4py import MPI


class Particle:
    def __init__(self, x0):
        self.position_i = []          # particle position
        self.velocity_i = []          # particle velocity
        self.pos_best_i = []          # best position individual
        self.err_best_i = -1          # best error individual
        self.err_i = -1               # error individual

        for i in range(0, num_dimensions):
            self.velocity_i.append(random.uniform(-1, 1))
            self.position_i.append(x0[i])

    # evaluate current fitness
    def evaluate(self, costFunc):
        self.err_i = costFunc(self.position_i)

        # check to see if the current position is an individual best
        if self.err_i < self.err_best_i or self.err_best_i == -1:
            self.pos_best_i = self.position_i.copy()
            self.err_best_i = self.err_i

    # update new particle velocity
    def update_velocity(self, pos_best_g):
        # constant inertia weight (how much to weigh the previous velocity)
        w = 0.5
        c1 = 1        # cognative constant
        c2 = 2        # social constant

        for i in range(0, num_dimensions):
            r1 = random.random()
            r2 = random.random()

            vel_cognitive = c1*r1*(self.pos_best_i[i]-self.position_i[i])
            vel_social = c2*r2*(pos_best_g[i]-self.position_i[i])
            self.velocity_i[i] = w*self.velocity_i[i]+vel_cognitive+vel_social

    # update the particle position based off new velocity updates
    def update_position(self, bounds):
        for i in range(0, num_dimensions):
            self.position_i[i] = self.position_i[i]+self.velocity_i[i]

            # adjust maximum position if necessary
            if self.position_i[i] > bounds[i][1]:
                self.position_i[i] = bounds[i][1]

            # adjust minimum position if neseccary
            if self.position_i[i] < bounds[i][0]:
                self.position_i[i] = bounds[i][0]


class PSO():
    def __init__(self, costFunc, x0, bounds, num_particles, maxiter, verbose=False):
        global num_dimensions

        num_dimensions = len(x0[0])
        err_best_g = -1                   # best error for group
        pos_best_g = []                   # best position for group

        # establish the swarm
        swarm = []
        for i in range(0, num_particles):
            swarm.append(Particle(x0[i]))

        # begin optimization loop
        i = 0
        while i < maxiter:
            if verbose:
                print(f'iter: {i:>4d}, best solution: {err_best_g:10.6f}')
            # cycle through particles in swarm and evaluate fitness
            for j in range(0, num_particles):
                swarm[j].evaluate(costFunc)

                # determine if current particle is the best (globally)
                if swarm[j].err_i < err_best_g or err_best_g == -1:
                    pos_best_g = list(swarm[j].position_i)
                    err_best_g = float(swarm[j].err_i)
                
            # cycle through swarm and update velocities and position
            for j in range(0, num_particles):
                swarm[j].update_velocity(pos_best_g)
                swarm[j].update_position(bounds)
            i += 1
            print(pos_best_g)
            print(err_best_g)

        self.err_best_g=err_best_g
        self.pos_best_g=pos_best_g
        # print final results
        print('\nFINAL SOLUTION:')
        print(f'   > {pos_best_g}')
        print(f'   > {err_best_g}\n')

# #--- RUN ----------------------------------------------------------------------+


# initial = [5, 5]               # initial starting location [x1,x2...]
# # input bounds [(x1_min,x1_max),(x2_min,x2_max)...]
# bounds = [(-10, 10), (-10, 10)]
# PSO(func1, initial, bounds, num_particles=15, maxiter=30, verbose=True)

# #--- END ----------------------------------------------------------------------+

def run_circuit(qureg, theta):
    """ 
    Runs target circuit using the theta.
    Args:
        qureg (allocate_qureg):Assigned qubit register
        theta: Parameters that need to be optimized.
    Returns:
        None.
    """
    X | qureg[0]
    X | qureg[1]
    X | qureg[2]
    X | qureg[3]
    X | qureg[4]
    X | qureg[5]
    X | qureg[6]
    X | qureg[7]
    X | qureg[8]
    X | qureg[9]
    Rx(theta[0]) | qureg[0]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[0]
    H | qureg[6]
    H | qureg[7]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    H | qureg[7]
    H | qureg[6]
    Rx(theta[0]) | qureg[3]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[0]
    H | qureg[10]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    Rz(theta[0]) | qureg[10]
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[10]
    Rx(theta[2]) | qureg[0]
    H | qureg[0]
    H | qureg[1]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    H | qureg[1]
    H | qureg[0]
    Rx(theta[0]) | qureg[8]
    Rx(theta[0]) | qureg[9]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[9]
    Rx(theta[2]) | qureg[8]
    Rx(theta[0]) | qureg[2]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[8]
    H | qureg[9]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[9]
    H | qureg[8]
    H | qureg[2]
    H | qureg[3]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[8]
    H | qureg[9]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[13]
    H | qureg[12]
    H | qureg[9]
    Rx(theta[2]) | qureg[8]
    H | qureg[2]
    Rx(theta[0]) | qureg[3]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[1]
    H | qureg[2]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    H | qureg[1]
    H | qureg[4]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    H | qureg[4]
    H | qureg[1]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[6]
    H | qureg[7]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[13]
    H | qureg[12]
    H | qureg[7]
    Rx(theta[2]) | qureg[6]
    H | qureg[1]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[1]
    Rx(theta[0]) | qureg[2]
    H | qureg[3]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[11]
    H | qureg[10]
    H | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[8]
    H | qureg[9]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[9]
    Rx(theta[2]) | qureg[8]
    H | qureg[0]
    H | qureg[1]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    H | qureg[1]
    H | qureg[0]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[1]
    Rx(theta[2]) | qureg[0]
    H | qureg[2]
    H | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    H | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[2]
    H | qureg[1]
    H | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    H | qureg[11]
    H | qureg[5]
    H | qureg[1]
    Rx(theta[0]) | qureg[1]
    H | qureg[4]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    H | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[1]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[0]
    H | qureg[4]
    H | qureg[10]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    H | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[0]
    H | qureg[8]
    Rx(theta[0]) | qureg[9]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[9]
    H | qureg[8]
    H | qureg[4]
    H | qureg[5]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    H | qureg[5]
    H | qureg[4]
    H | qureg[2]
    H | qureg[3]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[3]
    H | qureg[2]
    H | qureg[6]
    Rx(theta[0]) | qureg[7]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[7]
    H | qureg[6]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[0]
    H | qureg[1]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    H | qureg[12]
    H | qureg[1]
    Rx(theta[2]) | qureg[0]
    H | qureg[2]
    H | qureg[3]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    H | qureg[3]
    H | qureg[2]
    H | qureg[3]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[3]
    H | qureg[1]
    H | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    H | qureg[4]
    H | qureg[1]
    H | qureg[2]
    H | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    H | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[2]
    H | qureg[5]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    H | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[2]
    H | qureg[3]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    H | qureg[12]
    H | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[4]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[4]
    H | qureg[3]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    H | qureg[8]
    H | qureg[9]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[1]
    Rx(theta[2]) | qureg[0]
    H | qureg[3]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[3]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[1]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    H | qureg[3]
    H | qureg[2]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    H | qureg[2]
    Rx(theta[0]) | qureg[1]
    H | qureg[4]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[12]
    H | qureg[11]
    H | qureg[4]
    Rx(theta[2]) | qureg[1]
    H | qureg[3]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[11]
    H | qureg[3]
    H | qureg[4]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[12]
    H | qureg[4]
    H | qureg[4]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[5]
    H | qureg[4]
    Rx(theta[0]) | qureg[6]
    Rx(theta[0]) | qureg[7]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[7]
    Rx(theta[2]) | qureg[6]
    H | qureg[2]
    H | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    H | qureg[10]
    H | qureg[4]
    H | qureg[2]
    H | qureg[0]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[1]
    H | qureg[0]
    Rx(theta[0]) | qureg[6]
    H | qureg[7]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[7]
    Rx(theta[2]) | qureg[6]
    Rx(theta[0]) | qureg[0]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[0]
    H | qureg[8]
    Rx(theta[0]) | qureg[9]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[6]
    H | qureg[7]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[7]
    Rx(theta[2]) | qureg[6]
    H | qureg[0]
    Rx(theta[0]) | qureg[10]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    Rz(theta[1]) | qureg[10]
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[10]
    H | qureg[0]
    H | qureg[0]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[1]
    H | qureg[0]
    H | qureg[3]
    H | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    H | qureg[11]
    H | qureg[5]
    H | qureg[3]
    Rx(theta[0]) | qureg[2]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[2]
    H | qureg[8]
    Rx(theta[0]) | qureg[9]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[1]
    H | qureg[3]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    H | qureg[3]
    Rx(theta[0]) | qureg[6]
    Rx(theta[0]) | qureg[7]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[7]
    Rx(theta[2]) | qureg[6]
    Rx(theta[0]) | qureg[1]
    H | qureg[4]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[4]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[4]
    H | qureg[5]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[13]
    H | qureg[12]
    H | qureg[5]
    Rx(theta[2]) | qureg[4]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[0]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[3]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[0]
    Rx(theta[0]) | qureg[4]
    H | qureg[5]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[5]
    Rx(theta[2]) | qureg[4]
    Rx(theta[0]) | qureg[3]
    H | qureg[4]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    H | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[2]
    H | qureg[0]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[3]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    H | qureg[0]
    H | qureg[2]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[2]
    H | qureg[6]
    Rx(theta[0]) | qureg[7]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[7]
    H | qureg[6]
    Rx(theta[0]) | qureg[0]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[3]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[1]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[1]
    H | qureg[3]
    H | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    H | qureg[4]
    H | qureg[3]
    Rx(theta[0]) | qureg[6]
    Rx(theta[0]) | qureg[7]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[7]
    Rx(theta[2]) | qureg[6]
    Rx(theta[0]) | qureg[1]
    H | qureg[5]
    H | qureg[11]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    H | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[3]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[0]
    H | qureg[1]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[11]
    H | qureg[10]
    H | qureg[1]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[0]
    H | qureg[5]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    H | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[1]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[1]
    Rx(theta[2]) | qureg[0]
    H | qureg[5]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    Rx(theta[2]) | qureg[13]
    H | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[3]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[6]
    H | qureg[7]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[11]
    H | qureg[10]
    H | qureg[7]
    Rx(theta[2]) | qureg[6]
    Rx(theta[0]) | qureg[4]
    H | qureg[12]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[12]
    Rx(theta[2]) | qureg[4]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[4]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[3]
    H | qureg[5]
    H | qureg[11]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    H | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[0]
    H | qureg[2]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[2]
    H | qureg[10]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    Rz(theta[3]) | qureg[10]
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[10]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[1]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[8]
    H | qureg[9]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[9]
    Rx(theta[2]) | qureg[8]
    H | qureg[2]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    H | qureg[2]
    H | qureg[8]
    Rx(theta[0]) | qureg[9]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[3]
    H | qureg[11]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[3]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[11]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[12]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[2]
    H | qureg[4]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[4]
    Rx(theta[0]) | qureg[2]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[2]
    Rx(theta[0]) | qureg[10]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    Rz(theta[2]) | qureg[10]
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[10]
    H | qureg[2]
    Rx(theta[0]) | qureg[4]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[2]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[0]
    H | qureg[0]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    H | qureg[0]
    H | qureg[3]
    H | qureg[4]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    H | qureg[4]
    H | qureg[3]
    H | qureg[0]
    H | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    H | qureg[5]
    H | qureg[0]
    Rx(theta[0]) | qureg[2]
    H | qureg[5]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    H | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[6]
    Rx(theta[0]) | qureg[7]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[7]
    H | qureg[6]
    Rx(theta[0]) | qureg[8]
    Rx(theta[0]) | qureg[9]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[9]
    Rx(theta[2]) | qureg[8]
    Rx(theta[0]) | qureg[0]
    H | qureg[5]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    H | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[3]
    H | qureg[4]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[4]
    Rx(theta[2]) | qureg[3]
    H | qureg[4]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    H | qureg[4]
    H | qureg[1]
    H | qureg[4]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[4]
    H | qureg[1]
    H | qureg[4]
    Rx(theta[0]) | qureg[5]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[5]
    H | qureg[4]
    H | qureg[0]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[0]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[4]
    H | qureg[1]
    H | qureg[2]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[3]
    H | qureg[2]
    Rx(theta[0]) | qureg[2]
    H | qureg[3]
    Rx(theta[0]) | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[1]
    H | qureg[11]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[11]
    Rx(theta[2]) | qureg[1]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[3]
    H | qureg[4]
    H | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    H | qureg[5]
    H | qureg[4]
    Rx(theta[0]) | qureg[4]
    H | qureg[5]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[11]
    H | qureg[10]
    H | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[0]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[0]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[6]
    H | qureg[7]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    H | qureg[7]
    H | qureg[6]
    H | qureg[6]
    H | qureg[7]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[7]
    H | qureg[6]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[5]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    Rx(theta[0]) | qureg[4]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    H | qureg[6]
    Rx(theta[0]) | qureg[7]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[7]
    H | qureg[6]
    H | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[1]
    H | qureg[0]
    H | qureg[2]
    H | qureg[3]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    H | qureg[3]
    H | qureg[2]
    H | qureg[1]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[1]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[0]
    H | qureg[0]
    Rx(theta[0]) | qureg[1]
    H | qureg[12]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[1], qureg[12] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[1]
    H | qureg[0]
    H | qureg[0]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    H | qureg[0]
    H | qureg[4]
    H | qureg[5]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[5], qureg[12] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    H | qureg[5]
    H | qureg[4]
    Rx(theta[0]) | qureg[1]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    Rx(theta[2]) | qureg[1]
    H | qureg[3]
    H | qureg[4]
    Rx(theta[0]) | qureg[11]
    H | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[12]
    Rx(theta[2]) | qureg[11]
    H | qureg[4]
    H | qureg[3]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[3]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[3]
    Rx(theta[2]) | qureg[2]
    H | qureg[1]
    Rx(theta[0]) | qureg[4]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[4]
    H | qureg[1]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[3]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[0]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[3], qureg[12] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[3]
    Rx(theta[2]) | qureg[2]
    Rx(theta[0]) | qureg[6]
    Rx(theta[0]) | qureg[7]
    Rx(theta[0]) | qureg[12]
    H | qureg[13]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[7], qureg[12] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[13]
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[7]
    Rx(theta[2]) | qureg[6]
    H | qureg[0]
    Rx(theta[0]) | qureg[5]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[2]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[5]
    H | qureg[0]
    H | qureg[3]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    H | qureg[0]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[0]
    H | qureg[8]
    H | qureg[9]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    H | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[4]
    H | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[4]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[1]
    H | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[11]
    H | qureg[1]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[4]
    H | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[12]
    H | qureg[11]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[3]
    H | qureg[1]
    Rx(theta[0]) | qureg[5]
    H | qureg[11]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    H | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[1]
    H | qureg[8]
    H | qureg[9]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    H | qureg[9]
    H | qureg[8]
    Rx(theta[0]) | qureg[3]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    Rx(theta[2]) | qureg[3]
    Rx(theta[0]) | qureg[8]
    Rx(theta[0]) | qureg[9]
    H | qureg[12]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[9], qureg[12] )
    CX | ( qureg[8], qureg[9] )
    Rx(theta[2]) | qureg[13]
    H | qureg[12]
    Rx(theta[2]) | qureg[9]
    Rx(theta[2]) | qureg[8]
    Rx(theta[0]) | qureg[0]
    Rx(theta[0]) | qureg[1]
    H | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[1]
    Rx(theta[2]) | qureg[0]
    Rx(theta[0]) | qureg[2]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[5]
    Rx(theta[2]) | qureg[2]
    H | qureg[2]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[2]
    Rx(theta[0]) | qureg[5]
    H | qureg[13]
    CX | ( qureg[5], qureg[6] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[7], qureg[8] )
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[5], qureg[6] )
    H | qureg[13]
    Rx(theta[2]) | qureg[5]
    Rx(theta[0]) | qureg[8]
    Rx(theta[0]) | qureg[9]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[9]
    Rx(theta[2]) | qureg[8]
    H | qureg[0]
    H | qureg[4]
    H | qureg[10]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[3]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    H | qureg[10]
    H | qureg[4]
    H | qureg[0]
    H | qureg[0]
    H | qureg[5]
    H | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[13]
    H | qureg[10]
    H | qureg[5]
    H | qureg[0]
    Rx(theta[0]) | qureg[3]
    H | qureg[4]
    H | qureg[11]
    H | qureg[12]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[4]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[4], qureg[11] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[12]
    H | qureg[11]
    H | qureg[4]
    Rx(theta[2]) | qureg[3]
    H | qureg[3]
    Rx(theta[0]) | qureg[4]
    H | qureg[10]
    H | qureg[13]
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    H | qureg[13]
    H | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[3]
    H | qureg[1]
    H | qureg[4]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[4]
    H | qureg[1]
    H | qureg[1]
    H | qureg[5]
    Rx(theta[0]) | qureg[11]
    H | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[3]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    H | qureg[13]
    Rx(theta[2]) | qureg[11]
    H | qureg[5]
    H | qureg[1]
    H | qureg[0]
    Rx(theta[0]) | qureg[5]
    Rx(theta[0]) | qureg[11]
    Rx(theta[0]) | qureg[12]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[1]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[5], qureg[11] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[0], qureg[1] )
    Rx(theta[2]) | qureg[12]
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[5]
    H | qureg[0]
    Rx(theta[0]) | qureg[2]
    H | qureg[4]
    H | qureg[10]
    H | qureg[12]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    Rz(theta[2]) | qureg[12]
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[12]
    H | qureg[10]
    H | qureg[4]
    Rx(theta[2]) | qureg[2]
    H | qureg[2]
    H | qureg[5]
    Rx(theta[0]) | qureg[10]
    H | qureg[13]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[4]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[5], qureg[10] )
    CX | ( qureg[4], qureg[5] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    H | qureg[13]
    Rx(theta[2]) | qureg[10]
    H | qureg[5]
    H | qureg[2]
    Rx(theta[0]) | qureg[8]
    H | qureg[9]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[8], qureg[9] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[0]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[9], qureg[10] )
    CX | ( qureg[8], qureg[9] )
    H | qureg[11]
    H | qureg[10]
    H | qureg[9]
    Rx(theta[2]) | qureg[8]
    Rx(theta[0]) | qureg[2]
    H | qureg[3]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[11]
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[3], qureg[10] )
    CX | ( qureg[2], qureg[3] )
    Rx(theta[2]) | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[3]
    Rx(theta[2]) | qureg[2]
    H | qureg[6]
    H | qureg[7]
    Rx(theta[0]) | qureg[10]
    H | qureg[11]
    CX | ( qureg[6], qureg[7] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[2]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[7], qureg[10] )
    CX | ( qureg[6], qureg[7] )
    H | qureg[11]
    Rx(theta[2]) | qureg[10]
    H | qureg[7]
    H | qureg[6]
    H | qureg[0]
    Rx(theta[0]) | qureg[1]
    H | qureg[10]
    H | qureg[11]
    CX | ( qureg[0], qureg[1] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    Rz(theta[1]) | qureg[11]
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[1], qureg[10] )
    CX | ( qureg[0], qureg[1] )
    H | qureg[11]
    H | qureg[10]
    Rx(theta[2]) | qureg[1]
    H | qureg[0]
    H | qureg[1]
    Rx(theta[0]) | qureg[4]
    Rx(theta[0]) | qureg[10]
    Rx(theta[0]) | qureg[13]
    CX | ( qureg[1], qureg[2] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[12], qureg[13] )
    Rz(theta[1]) | qureg[13]
    CX | ( qureg[12], qureg[13] )
    CX | ( qureg[11], qureg[12] )
    CX | ( qureg[10], qureg[11] )
    CX | ( qureg[4], qureg[10] )
    CX | ( qureg[3], qureg[4] )
    CX | ( qureg[2], qureg[3] )
    CX | ( qureg[1], qureg[2] )
    Rx(theta[2]) | qureg[13]
    Rx(theta[2]) | qureg[10]
    Rx(theta[2]) | qureg[4]
    H | qureg[1]

def cost(theta, target_state, eng):
    n=len(target_state)
    qureg = eng.allocate_qureg(n)
    run_circuit(qureg,theta)
    
    eng.flush()
    count=eng.backend.get_probability(target_state,qureg)
    All(Measure) | qureg

    return 1-count


def pre_cal_param(eng,final_state,N=20):
    opt_cost=1
    opt_param=[]
    for i in range(N):
        rand_theta=[np.pi*2*np.random.rand() for i in range(5)]
        rand_cost=cost(rand_theta,final_state,eng)
        if rand_cost<opt_cost:
            opt_cost=rand_cost
            opt_param=rand_theta
    return opt_param




def calculate_theta(eng, final_state):
    """
    Calculate the optimal theta value.
    Args:
    eng (MainEngine): Main compiler engine to use.
    final_state (list): list of quantum bit states.
    Returns:
    list of best theta values.
    """
    n_particle=30
#    theta0=pre_cal_param(eng,final_state)
    theta0=[[np.pi*2*np.random.rand() for i in range(5)] for j in range(n_particle)]
    boundary=[(0,np.pi*2) for t in theta0]
    # result=minimize(cost,x0=theta,args=(final_state,eng),method="COBYLA")
    def warpped_cost(theta):
        return cost(theta, final_state, eng)

    result=PSO(warpped_cost,theta0,boundary,n_particle,maxiter=50,verbose=True)
    
    return result.pos_best_g


if __name__ == "__main__":
    import datetime
    
    backend = SimulatorMPI(gate_fusion=True,num_local_qubits=14)
    cache_depth = 2000
    rule_set = DecompositionRuleSet(modules=[rules])
    engines = [TagRemover()
    , LocalOptimizer(cache_depth)
    , AutoReplacer(rule_set)
    , TagRemover()
    , LocalOptimizer(cache_depth)
    , GreedyScheduler()
    ]

    eng = HiQMainEngine(backend, engines)
    t1=datetime.datetime.now()
    print(t1)
    
    final_state = [0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1]

    theta = calculate_theta(eng, final_state)
    # theta=np.ones(5)
    
    qureg = eng.allocate_qureg(14)
    run_circuit(qureg, theta)
    eng.flush()

    t2=datetime.datetime.now()
    print(t2)
    print(t2-t1)
    
    print("probability")
    print(eng.backend.get_probability(final_state, qureg))
    All(Measure) | qureg