#
#  Copyright 2022 Jorge Velazquez Castro
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
from epinetwork.Models import VectorBorne
from epinetwork.MobilityNetwork import MobilityNetwork
from epinetwork.Control_protocol import controlProtocol
from epinetwork.SimulationsEnsemble import simulationsEnsemble
from epinetwork.Control_protocol import RandomControl, IndexBasedControl, noControl

import matplotlib
import matplotlib.pyplot as plt

def detailed():
    #parameteres
    number_of_simulations = 100
    n = 6 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 0.67
    param[1] = gamma = 1./7.
    param[2] = beta_v = 0.67  #5
    param[3] = mu_v = 1./8.
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 1500.
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 200.
    W = y[4] = 0.#

    P = MobilityNetwork()
    #P.barabsi_albert(n,m=2,min_residential=min_residential)
    P.binomial(n,b,min_residential=min_residential)
    vectorModel = VectorBorne(n,params=param, network=P)

    sim = simulation(vectorModel)
    sim.set_initial_conditions_all_patches(y)
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(400) #How many "days" to simulate

    sim.model.set_protocol(IndexBasedControl(update_interval=None))
    index_sim = sim.copy()
    index_sim.run()

    sim.model.set_protocol(RandomControl(update_interval = None))
    rnd_sim = sim.copy()
    rnd_sim.run()

    sim.model.set_protocol(noControl)
    no_control_sim = sim.run().copy()

    fig, ax1 = plt.subplots()
    t = index_sim.time
    for i in range(index_sim.No_patches):
        ax1.plot(t,index_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('Index control')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t = rnd_sim.time
    for i in range(rnd_sim.No_patches):
        ax1.plot(t,rnd_sim.evolution[i,1,:],label="Patch " + str(i))
    ax1.set_title('Random control')
    ax1.legend()

    fig, ax1 = plt.subplots()
    t = no_control_sim.time
    for i in range(no_control_sim.No_patches):
        ax1.plot(t,no_control_sim.evolution[i,1,:],label="Patch " + str(i))
        print("Patch " + str(i) + ":")
        print("Final size " + str(no_control_sim.evolution[i,2,-1]))
        print("Index " + str(index_sim.model.control.indices[i]))
    ax1.set_title('No control')
    ax1.legend()

    fig, ax1 = plt.subplots()
    ax1.plot(t,no_control_sim.evolution_totals[1,:], label="No control")
    ax1.plot(t,rnd_sim.evolution_totals[1,:], label="Random control")
    ax1.plot(t,index_sim.evolution_totals[1,:], label="Index control")
    ax1.set_title('Total Infected')
    ax1.legend()


    plt.show()

def index():
    n = 1 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 1./10. # 0.67
    param[1] = gamma = 1./7. #1./7.
    param[2] = beta_v = 1./10. #0.67  #5
    param[3] = mu_v = 1./10. # 1./8.
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 15000. #15000
    W = y[4] = 0.#

    #P = MobilityNetwork()
    #P.barabsi_albert(n,m=2,min_residential=min_residential)
    #P.binomial(n,b,min_residential=min_residential)
    vectorModel = VectorBorne(n,params=param)
    #vectorModel.p.matrix[0,1]=0.1
    #vectorModel.p.matrix[1,0]=0.0
    #vectorModel.p.matrix[1,1]=1.0
    #vectorModel.p.matrix[0,0]=0.9


    sim = simulation(vectorModel)
    sim.set_initial_conditions_all_patches(y)
    N = sim.node[:,0]
    Nv = sim.node[:,3]
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(0,y) #Se establece esta condicion inicial en la zona 1
     #How many "days" to simulate

    R_index_list = []
    R_sim_list = []

    beta_range = np.arange(0.1,6.,0.1) #0.01
    for betav in beta_range :
        param[2]=betav
        #param[0]=betav
        vectorModel.set_patches_params(param)
        sim.set_simulation_time(int(1000000 - 1000*betav))
        R_index_list.append(vectorModel.local_final_size(N,N,Nv).copy())
        sim.run()
        R_sim_list.append(sim.evolution[:,2,-1].copy())

    R_index = np.array(R_index_list)
    R_sim = np.array(R_sim_list)
    for patch in range(vectorModel.number_of_patches):
         fig, ax1 = plt.subplots()
         #ax1.plot(beta_range,R_index[:,patch],label=r'$R_{approx}$')
         #ax1.plot(beta_range,R_sim[:,patch], label=r'$R_{exact}$')
         ax1.plot(mu_v/beta_range,R_index[:,patch],label=r'$R_{approx}$')
         ax1.plot(mu_v/beta_range,R_sim[:,patch], label=r'$R_{exact}$')
         ax1.legend()
         ax1.set_xlabel(r'$\mu/\beta_v$')
         ax1.set_ylabel(r'$R(\infty)$')
         #ax1.set_title(r'$R(\infty)$ vs  $\mu/\beta$' )
         fig, ax2 = plt.subplots()
         R_sim_vec=R_sim[:,patch].flatten()
         R_index_vec=R_index[:,patch].flatten()
         Error =[ (R_in_index - R_in_sim)/R_in_sim if R_in_index >1. else 0. for R_in_index,R_in_sim in zip(R_index_vec,R_sim_vec)]
         ax2.plot(mu_v/beta_range,Error)
         #ax1.plot(mu_v/beta_range,R_sim[:,patch], label=r'$R_{exact}$')
         ax2.legend()
         ax2.set_xlabel(r'$\mu/\beta_v$')
         ax2.set_ylabel(r'$|\frac{R_{approx} - R_{exact}}{R_{exact}}|$')


    plt.show()


def main():
    #parameteres
    number_of_simulations = 50
    n = 6 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 0.67
    param[1] = gamma = 1./7.
    param[2] = beta_v =5
    param[3] = mu_v = 1./8.
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 1500.
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 200.
    W = y[4] = 0.#

    P = MobilityNetwork()
    P.binomial(n,b,min_residential)
    control = IndexBasedControl(update_interval=None)
    vectorModel = VectorBorne(params=param, network=P, control=control)

    sim = simulation(vectorModel)
    sim.set_initial_conditions_all_patches(y)
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(80) #How many "days" to simulate

    ensemble = simulationsEnsemble()
    ensemble.add_simulation(sim)

    for i in range(number_of_simulations):
        P.binomial(n,b,min_residential)
        ensemble.add_simulation(sim)

    ensemble.run_all_simulations()

    ensemble.plot_infected_average()
    ensemble.plot_recovered_average()
    plt.show()

def comparison():
    #parameteres
    number_of_simulations = 100
    n = 6 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 0.67
    param[1] = gamma = 1./7.
    param[2] = beta_v = 0.67  #5
    param[3] = mu_v = 1./8.
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 1500.
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 200.
    W = y[4] = 0.#

    P = MobilityNetwork()
    vectorModel = VectorBorne(n,params=param, network=P)

    sim = simulation(vectorModel)
    sim.set_initial_conditions_all_patches(y)
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
    sim.set_simulation_time(400) #How many "days" to simulate

    index_control = IndexBasedControl(update_interval=None)
    rand_control = RandomControl(update_interval = None)
    no_control = noControl()

    IndexControledEnsemble = simulationsEnsemble()
    RandControledEnsemble = simulationsEnsemble()
    NoControlEnsemble = simulationsEnsemble()

    for i in range(number_of_simulations):
        #Both ensembles have the same Mobility Network P
        #P.barabsi_albert(n,m=2,min_residential=min_residential)
        sim.model.p.barabsi_albert(n,m=2,min_residential=min_residential)
        #otros tipos de redes:
        #P.binomial(n,b,min_residential)
        #P.newman_watts_strogatz(n,k,p,min_residential)
        sim.model.set_protocol(index_control)
        IndexControledEnsemble.add_simulation(sim) #when a simulation is addet to an ensemble it makes a copy, thus we car resuse the object simulation
        sim.model.set_protocol(rand_control)
        RandControledEnsemble.add_simulation(sim)
        sim.model.set_protocol(no_control)
        NoControlEnsemble.add_simulation(sim)

    IndexControledEnsemble.multiprocessing_run_all(No_of_process=8)
    RandControledEnsemble.multiprocessing_run_all(No_of_process=8)
    NoControlEnsemble.multiprocessing_run_all(No_of_process=8)

    y_index_controlled, t = IndexControledEnsemble.get_evolution_average()
    y_rnd_controlled, t = RandControledEnsemble.get_evolution_average()
    y, t = NoControlEnsemble.get_evolution_average()

    fig, ax1 = plt.subplots()
    ax1.plot(t,y_index_controlled[1],label="Index controlled")
    ax1.plot(t,y_rnd_controlled[1],label="Random controlled")
    ax1.plot(t,y[1],label="Not controlled")
    ax1.set_ylabel('Infected population')
    ax1.legend()

    fig2, ax2 = plt.subplots()
    ax2.plot(t,y_index_controlled[4],label="Index controlled")
    ax2.plot(t,y_rnd_controlled[4],label="Random controlled")
    ax2.plot(t,y[4],label="Not controlled")
    ax2.set_ylabel('Invected vector population')
    ax2.legend()

    plt.show()

if (True):
    #main()
    #comparison()
    #detailed()
    index()
