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

def index(plot=True):
    n = 4 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 0.0686 # (0.0146,0.2241)
    param[1] = gamma = 1./7. #1./7.
    param[2] = beta_v = 0.4307 # (0.1299, 1.6821)
    param[3] = mu_v = 1./15. #
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 15000. #15000
    W = y[4] = 0.#

    P = MobilityNetwork()
    P.barabsi_albert(n,m=2,min_residential=min_residential)
    #P.binomial(n,b,min_residential=min_residential)
    vectorModel = VectorBorne(n,params=param,network=P)
    rng = np.random.default_rng()
    vectorModel.beta_v = (1.68-0.13)*rng.random((n)) + 0.13
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

    R_index_inf_list = []
    R_index_sup_list = []
    R_sim_list = []

    beta_range = np.arange(0.13,1.6,0.01)
    for betav in beta_range :
        print("calculando " , betav)
        param[2]=betav
        #param[0]=betav
        vectorModel.set_patches_params(param)
        sim.set_simulation_time(int(10000000))
        sim.set_early_stopping()
        R_index_inf_list.append(vectorModel.local_final_size_inf(N,N,Nv).copy()/N)
        R_index_sup_list.append(vectorModel.local_final_size_sup(N,N,Nv).copy()/N)
        sim.run()
        R_sim_list.append(sim.evolution[:,2,-1].copy()/N)

    R_inf = np.array(R_index_inf_list)
    R_sup = np.array(R_index_sup_list)
    R_sim = np.array(R_sim_list)
    if(plot==True):
        fig, axs = plt.subplots(2,2)
        ax = axs.flat
        for patch in range(vectorModel.number_of_patches):
            ax[patch].plot(beta_range,R_inf[:,patch],label=r'$R_{sup}$',marker='*')
            ax[patch].plot(beta_range,R_sup[:,patch],label=r'$R_{inf}$',marker='x')
            ax[patch].plot(beta_range,R_sim[:,patch], label=r'$R_{exact}$')
            #ax1.plot(mu_v/beta_range,R_inf[:,patch],label=r'$R_{inf}$')
            #ax1.plot(mu_v/beta_range,R_inf[:,patch],label=r'$R_{inf}$')
            #ax1.plot(mu_v/beta_range,R_sim[:,patch], label=r'$R_{exact}$')
            ax[patch].legend()
            ax[patch].set_xlabel(r'$\beta_v$')
            ax[patch].set_ylabel(r'$R(\infty)$')
            ax[patch].label_outer()
            #ax1.set_title(r'$R(\infty)$ vs  $\mu/\beta$' )
            #fig, ax2 = plt.subplots()
            #R_sim_vec=R_sim[:,patch].flatten()
            #R_index_vec=R_index[:,patch].flatten()
            #Error =[ (R_in_index - R_in_sim)/R_in_sim if R_in_index >1. else 0. for R_in_index,R_in_sim in zip(R_index_vec,R_sim_vec)]
            #ax2.plot(mu_v/beta_range,Error)
            #ax1.plot(mu_v/beta_range,R_sim[:,patch], label=r'$R_{exact}$')
            #ax2.legend()
            #ax2.set_xlabel(r'$\mu/\beta_v$')
            #ax2.set_ylabel(r'$|\frac{R_{approx} - R_{exact}}{R_{exact}}|$')
        plt.show()
    
    return np.array([R_sim,R_inf,R_sup]) , beta_range

def index_error():
    R_s , beta = index(plot=False)
    Err_sup = np.max(np.abs(R_s[1,:,:]-R_s[0,:,:]),axis=-1)
    Err_inf = np.max(np.abs(R_s[2,:,:]-R_s[0,:,:]),axis=-1)
    plt.plot(beta,Err_sup,label=r'$R_{sup}$')
    plt.plot(beta,Err_inf,label=r'$R_{inf}$')
    plt.xlabel(r'$\beta_v$')
    plt.ylabel(r'$max|R_{approx} - R_{exact}|$')
    plt.legend()
    plt.show()


def index_order(plot=True):
    n = 25 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.6 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(5)
    param[0] = beta_h = 0.0686 # (0.0146,0.2241) #0.0686
    param[1] = gamma = 1./7. #1./7.
    param[2] = beta_v = 0.8 # (0.1299, 1.6821) sirve hasta 1.5 max para betha_h 0.068
    param[3] = mu_v = 1./15. #
    #initial conditions para una zona
    y = np.zeros(5)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0
    V = y[3] = 15000. #15000
    W = y[4] = 0.#

    P = MobilityNetwork()
    P.barabsi_albert(n,m=2,min_residential=min_residential)
    #P.binomial(n,b,min_residential=min_residential)
    vectorModel = VectorBorne(n,params=param,network=P)
    #rng = np.random.default_rng()
    #vectorModel.beta_v = (1.68-0.13)*rng.random((n)) + 0.13

    sim = simulation(vectorModel)
    #sim.set_initial_conditions_all_patches(y)
    rng = np.random.default_rng()
    S=rng.integers(low=15000,high=35000, size=n)

    for patch,s in enumerate(S):
        y[0] = s
        sim.set_initial_conditions_patch(patch,y)

    N = sim.node[:,0]
    Nv = sim.node[:,3]
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(0,y) #Se establece esta condicion inicial en la zona 1
     #How many "days" to simulate
    sim.set_simulation_time(2000000.)
    sim.set_early_stopping()
    R_index_inf = vectorModel.local_final_size_inf(N,N,Nv).copy()
    R_index_sup = vectorModel.local_final_size_sup(N,N,Nv).copy()
    sim.run()
    print(sim.t_events)
    R_sim = sim.y_events[:,2].copy()
    #ind_max = np.argmax(sim.y_events[:,1])

    #print(sim.t_events)
    #print(R_sim[ind_max])

    #sim.plot_all()
    #Normalization
    R_sim=R_sim/N
    R_index_inf=R_index_inf/N
    R_index_sup=R_index_sup/N


    R_sorted = np.array(sorted(zip(R_sim,R_index_inf,R_index_sup)))

    Err_sup=np.abs(R_sorted[:,1]-R_sorted[:,0])
    Err_inf=np.abs(R_sorted[:,2]-R_sorted[:,0])

    if plot==True:
        fig, ax = plt.subplots()
        x=range(n)
        ax.plot(x,R_sorted[:,0], label=r'$R_{exact}$')
        ax.scatter(x,R_sorted[:,1], label=r'$R_{inf}$',marker='*' )
        ax.scatter(x,R_sorted[:,2], label=r'$R_{sup}$',marker='x')
        ax.legend()
        ax.set_xlabel("zone")
        ax.set_ylabel(r'$R(\infty)$')
        ax.set_xticks(x)

        fig, ax2 = plt.subplots()
        ax2.plot(x,Err_sup,label=r'$R_{sup}$')
        ax2.plot(x,Err_inf,label=r'$R_{inf}$')
        ax2.legend()
        ax2.set_xlabel("zone")
        ax2.set_ylabel(r'$|R_{approx} - R_{exact}|$')
        ax2.set_xticks(x)

        plt.show()
    else:
        return R_sorted
    
def average_index_order():
    n=1000
    R = index_order(plot=False)
    R_sim=R[:,0]
    R_inf=R[:,1]
    R_sup=R[:,2]
    Err_sup = []
    Err_inf = []
    for sim in range(n):
        R = index_order(plot=False)
        R_sim+=R[:,0]
        R_inf+=R[:,1]
        R_sup+=R[:,2]
        Err_sup.append(R[:,1]-R[:,0])
        Err_inf.append(R[:,2]-R[:,0])
    R_sim=R_sim/float(n)
    R_inf=R_inf/float(n)
    R_sup=R_sup/float(n)
    Err_sup=np.array(Err_sup)
    Err_inf=np.array(Err_inf)

    #Normalization
    #R_sim=R_sim/np.max(R_sim)
    #R_inf=R_inf/np.max(R_inf)
    #R_sup=R_sup/np.max(R_sup)

    fig, ax = plt.subplots()
    x=range(len(R_sim))
    ax.scatter(x,R_sim, label=r'$R_{exact}$')
    ax.scatter(x,R_inf, label=r'$R_{sup}$',marker='*' )
    ax.scatter(x,R_sup, label=r'$R_{inf}$',marker='x')
    ax.legend()
    ax.set_xlabel("zone")
    ax.set_ylabel(r'$R(\infty)$')
    ax.set_xticks(x)

    fig, ax2 = plt.subplots()
    mean=np.mean(Err_sup,axis=0)
    lower_bar = np.percentile(Err_sup, 5, axis=0)
    upper_bar = np.percentile(Err_sup, 95, axis=0)
    #ax2.errorbar(x,mean,yerr=(lower_bar,upper_bar),lolims=lolims,uplims=uplims,linestyle='dotted')
    ax2.scatter(x,mean,label=r'$|R_{sup}-R_{exact}|$',marker='*',c='orange')
    ax2.fill_between(x, lower_bar, upper_bar, color='b', alpha=.1)
    ax2.set_xlabel("zone")
    #ax2.set_ylabel(r'$|R_{sup}-R_{exact}|$')
    ax2.set_xticks(x)

    #fig, ax3 = plt.subplots()
    mean=np.mean(Err_inf,axis=0)
    lower_bar = np.percentile(Err_inf, 5, axis=0)
    upper_bar = np.percentile(Err_inf, 95, axis=0)
    #ax3.errorbar(x,mean,yerr=(lower_bar,upper_bar),capsize=1.)
    ax2.scatter(x,-mean, label=r'$|R_{inf}-R_{exact}|$', marker='x',c='green')
    ax2.fill_between(x, -lower_bar, -upper_bar, color='b', alpha=.1)
    #ax3.set_xlabel("zone")
    #ax3.set_ylabel(r'$|R_{inf}-R_{exact}|$')
    #ax3.set_xticks(x)
    ax2.legend(loc='center right')
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
    #index()
    #index_error()
    index_order()
    #average_index_order()
