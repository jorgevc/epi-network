#
#  Copyright 2024 Jorge Velazquez Castro
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
import numpy as np

from epinetwork.Simulation import simulation
from epinetwork.Models import SEIR
from epinetwork.MobilityNetwork import MobilityNetwork
from epinetwork.SimulationsEnsemble import simulationsEnsemble
from epinetwork.Control_protocol import RandomControl, IndexBasedControl, noControl

import matplotlib
import matplotlib.pyplot as plt

def detailed():
    #parameteres
    number_of_simulations = 10
    n = 6 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(4)
    param[0] = beta = 0.67  # from 0.13 to 1.68
    param[1] = kappa = 1./7.
    param[2] = gamma = 0.67 
    param[3] = epsilon = 1./8.
    #initial conditions para una zona
    y = np.zeros(4)
    S = y[0] = 1500.  #from 1500 to 35000
    E = y[1] = 0.0
    I = y[2] = 0.0
    R = y[3] = 0.
    #cost able to spend
    cost = 6.

    P = MobilityNetwork()
    #P.barabsi_albert(n,m=2,min_residential=min_residential)
    P.binomial(n,b,min_residential=min_residential)
    seirModel = SEIR(n,params=param, network=P)
    
    ## Set random beta for each patch ##
    rng = np.random.default_rng()
    seirModel.beta = (1.68-0.13)*rng.random((n)) + 0.13

    sim = simulation(seirModel)
    ## Set random population in each patch ##
    rng = np.random.default_rng()
    S=rng.integers(low=15000,high=35000, size=n)
    for patch,s in enumerate(S):
        y[0] = s
        sim.set_initial_conditions_patch(patch,y)
    
    y[1]=y[1]+1. # Se agrega 1 expuesto a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(400) #How many "days" to simulate

    no_control_sim = sim.copy()
    no_control_sim.run()

    sim.model.set_protocol(IndexBasedControl(strategy='only_max',cost=cost,update_interval=None)) 
    only_max_sim = sim.copy()
    only_max_sim.run()

    sim.model.set_protocol(IndexBasedControl(strategy='equal_distributed',cost=cost,update_interval=None))
    equal_dist_sim = sim.copy()
    equal_dist_sim.run()

    sim.model.set_protocol(IndexBasedControl(strategy='cost_prop',cost=cost,update_interval=None))
    cost_prop_sim = sim.copy()
    cost_prop_sim.run()

    sim.model.set_protocol(IndexBasedControl(strategy='control_prop',cost=cost,update_interval=None))
    control_prop_sim = sim.copy()
    control_prop_sim.run()

    fig, ax1 = plt.subplots()
    t0 = no_control_sim.time
    for i in range(no_control_sim.No_patches):
        ax1.plot(t0,no_control_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('No control')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t1 = only_max_sim.time
    for i in range(only_max_sim.No_patches):
        ax1.plot(t1,only_max_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('Only max index controled')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t2 = equal_dist_sim.time
    for i in range(equal_dist_sim.No_patches):
        ax1.plot(t2,equal_dist_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('Equal distributed cost')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t3 = cost_prop_sim.time
    for i in range(cost_prop_sim.No_patches):
        ax1.plot(t3,cost_prop_sim.evolution[i,1,:],label="Patch " + str(i))
        print("Patch " + str(i) + ":")
        print("Final size " + str(cost_prop_sim.evolution[i,3,-1]))
        print("Index " + str(cost_prop_sim.model.control.indices[i]))
    ax1.set_title('Cost distributed proportional to index')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t4 = control_prop_sim.time
    for i in range(control_prop_sim.No_patches):
        ax1.plot(t4,control_prop_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('Control proportional to index')
    ax1.legend()

    fig, ax1 = plt.subplots()
    ax1.plot(t0,no_control_sim.evolution_totals[1,:], label="no control")
    ax1.plot(t1,only_max_sim.evolution_totals[1,:], label="only max index controled")
    ax1.plot(t2,equal_dist_sim.evolution_totals[1,:], label="equal distributed control")
    ax1.plot(t3,cost_prop_sim.evolution_totals[1,:], label="cost proportional to index")
    ax1.plot(t4,control_prop_sim.evolution_totals[1,:], label="control proportional to index")
    ax1.set_title('Total Infected')
    ax1.legend()


    plt.show()

if (True):
    detailed()