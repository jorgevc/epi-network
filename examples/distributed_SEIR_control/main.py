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
    param[0] = beta = 0.67  # Wuhan 0.59 to 1.68 ## GC 0.1 - 1.
    param[1] = kappa = 1./7.  #Wuhan 0.33        ## 0.07 - 0.5
    param[2] = gamma = 0.67  #Wuhan 0.2          ## 0.1 - 1
    param[3] = epsilon = 0.2 #porcentaje de infecciones asintomaticas
    #initial conditions para una zona
    y = np.zeros(4)
    S = y[0] = 1500.  #from 1500 to 35000
    E = y[1] = 0.0
    I = y[2] = 0.0
    R = y[3] = 0.
    #cost able to spend
    cost = 6.

    P = MobilityNetwork()
    P.barabsi_albert(n,m=2,min_residential=min_residential)
    #P.binomial(n,b,min_residential=min_residential)
    seirModel = SEIR(n,params=param, network=P)
    
    ## Set random beta for each patch ##
    rng = np.random.default_rng()
    #seirModel.beta = (1.68-0.13)*rng.random((n)) + 0.13
    seirModel.beta = (1.68-0.59)*rng.random((n)) + 0.59

    sim = simulation(seirModel)
    ## Set random population in each patch ##
    rng = np.random.default_rng()
    S=rng.integers(low=15000,high=35000, size=n)
    for patch,s in enumerate(S):
        y[0] = s
        sim.set_initial_conditions_patch(patch,y)
    
    y[1]=y[1]+1. # Se agrega 1 expuesto a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(250) #How many "days" to simulate

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
    ax1.set_title('Infected (Global count)')
    ax1.legend()


    plt.show()

def comparison():
    #parameteres
    number_of_simulations = 100
    n = 6 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(4)
    param[0] = beta = 0.67  # Wuhan 0.59 to 1.68 ## GC 0.1 - 1.
    param[1] = kappa = 1./7.  #Wuhan 0.33        ## 0.07 - 0.5
    param[2] = gamma = 0.67  #Wuhan 0.2          ## 0.1 - 1
    param[3] = epsilon = 0.2 #porcentaje de infecciones asintomaticas
    #initial conditions para una zona
    y = np.zeros(4)
    S = y[0] = 1500.  #from 1500 to 35000
    E = y[1] = 0.0
    I = y[2] = 0.0
    R = y[3] = 0.
    #cost able to spend
    cost = 6.

    #rng1 = np.random.default_rng()
    #rng2 = np.random.default_rng()

    P = MobilityNetwork()
    seirModel = SEIR(n,params=param, network=P)

    sim = simulation(seirModel)
    y[1]=y[1]+1. # Se agrega 1 expuesto a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(250) #How many "days" to simulate
    

    NoControlEnsemble = simulationsEnsemble()
    EqualDistEnsemble = simulationsEnsemble()
    CostPropEnsemble = simulationsEnsemble()
    CtrlPropEnsemble = simulationsEnsemble()
    for i in range(number_of_simulations):
        rng1 = np.random.default_rng(43+2*i)
        rng2 = np.random.default_rng(50121+3*i)
        ## Network ##
        #P.binomial(n,b,min_residential=min_residential)
        #sim.model.p.barabsi_albert(n,m=2,min_residential=min_residential)
        sim.model.p.binomial(n,b,min_residential=min_residential)

        ## Set random beta for each patch ##
        sim.model.beta = (1.68-0.59)*rng1.random((n)) + 0.59

        ## Set random population in each patch ##
        S=rng2.integers(low=15000,high=35000, size=n)
        for patch,s in enumerate(S):
            y[0] = s
            sim.set_initial_conditions_patch(patch,y)
    
        sim.model.set_protocol(noControl())
        NoControlEnsemble.add_simulation(sim)
        
        sim.model.set_protocol(IndexBasedControl(strategy='equal_distributed',cost=cost,update_interval=None))
        #sim.model.control.initial_time = 10.
        EqualDistEnsemble.add_simulation(sim)

        sim.model.set_protocol(IndexBasedControl(strategy='cost_prop',cost=cost,update_interval=None))
        #sim.model.control.initial_time = 10.
        CostPropEnsemble.add_simulation(sim)
   
        sim.model.set_protocol(IndexBasedControl(strategy='control_prop',cost=cost,update_interval=None))
        #sim.model.control.initial_time = 10.
        CtrlPropEnsemble.add_simulation(sim)
    
    NoControlEnsemble.multiprocessing_run_all(No_of_process=8)
    EqualDistEnsemble.multiprocessing_run_all(No_of_process=8)
    CostPropEnsemble.multiprocessing_run_all(No_of_process=8)
    CtrlPropEnsemble.multiprocessing_run_all(No_of_process=8)

    
    eq_cst_max=[]
    eq_cst_min=[]
    for eq_sim, cst_sim in zip(CtrlPropEnsemble.simulations,CostPropEnsemble.simulations):
        eq_maxs = np.max(eq_sim.evolution[:,2,:],axis=-1)
        eq_max = np.max(eq_maxs)
        eq_min = np.min(eq_maxs)
        cst_maxs = np.max(cst_sim.evolution[:,2,:],axis=-1)
        cst_max = np.max(cst_maxs)
        cst_min = np.min(cst_maxs)
        eq_cst_max.append(eq_max/cst_max)
        eq_cst_min.append(eq_min/cst_min)
    eq_cst_max = np.array(eq_cst_max)
    eq_cst_min = np.array(eq_cst_min)
    max_positive_rate = float(len(eq_cst_max[eq_cst_max>=1]))/float(len(eq_cst_max))
    min_positive_rate = float(len(eq_cst_min[eq_cst_min>=1]))/float(len(eq_cst_min))
    mean_max_variation = np.mean(eq_cst_max)
    mean_min_variation = np.mean(eq_cst_min)

    print("Percentage of sims that equal distributed is greater than cost dist:")
    print(max_positive_rate)
    print("Mean percentage of variation:")
    print(mean_max_variation)
    print("Percentage of sims that min of equal distributed is greater than cost dist:")
    print(min_positive_rate)
    print("Mean percentage of variation:")
    print(mean_min_variation)


    y_no_controlled, t = NoControlEnsemble.get_evolution_average()
    y_equal_controlled, t = EqualDistEnsemble.get_evolution_average()
    y_cost_prop, t = CostPropEnsemble.get_evolution_average()
    y_ctrl_prop, t = CtrlPropEnsemble.get_evolution_average()

    fig, ax1 = plt.subplots()
    ax1.plot(t,y_cost_prop[2],label="cost proportional to index")
    ax1.plot(t,y_ctrl_prop[2],label="control proportional to index")
    ax1.plot(t,y_equal_controlled[2],label="equal distributed control")
    ax1.plot(t,y_no_controlled[2],label="no control")
    ax1.set_ylabel('Infected (Ensemble average)')
    ax1.legend()

    fig2, ax2 = plt.subplots()
    ax2.plot(t,y_cost_prop[3],label="cost proportional to index")
    ax2.plot(t,y_ctrl_prop[3],label="control proportional to index")
    ax2.plot(t,y_equal_controlled[3],label="equal distributed control")
    ax2.plot(t,y_no_controlled[3],label="no control")
    ax2.set_ylabel('Removed (Cumulative Infections)')
    ax2.legend()

    plt.show()


if (True):
    #detailed()
    comparison()