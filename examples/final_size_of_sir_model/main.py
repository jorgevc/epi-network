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
from epinetwork.Models import SIR
from epinetwork.MobilityNetwork import MobilityNetwork
from epinetwork.Control_protocol import controlProtocol
from epinetwork.SimulationsEnsemble import simulationsEnsemble
from epinetwork.Control_protocol import RandomControl, IndexBasedControl, noControl

import matplotlib
import matplotlib.pyplot as plt

def index():
    n = 4 #numero de parches
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    param = np.zeros(2)
    param[0] = beta = 0.67 # 0.67
    param[1] = gamma = 1./7. #1./7.

    #initial conditions para una zona
    y = np.zeros(3)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0

    #P = MobilityNetwork()
    #P.barabsi_albert(n,m=2,min_residential=min_residential)
    #P.binomial(n,b,min_residential=min_residential)
    sirModel = SIR(n,params=param)
    P = MobilityNetwork()
    #P.barabsi_albert(n,m=2,min_residential=min_residential)
    b = 0.5 # parametro de la red binomial
    min_residential = 0.9 # diagonal de la matriz de mobilidad mayor a este numero
    P.binomial(n,b,min_residential=min_residential)
    sirModel.set_network(P)

    sim = simulation(sirModel)
    sim.set_initial_conditions_all_patches(y)
    N = sim.node[:,0]
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(0,y) #Se establece esta condicion inicial en la zona 1
    sim.run()
    sim.plot_all()

    R_inf = []
    iteration = range(1,10)

    for it in iteration:
        sirModel.final_size_max_iterations = it
        R_inf.append(sirModel.local_final_size(N,N).copy())
    R_inf = np.array(R_inf)

    fig, ax = plt.subplots(2,2)
    ax1=ax.flatten()
    for patch in range(sirModel.number_of_patches):

         ax1[patch].scatter(iteration,R_inf[:,patch])
         ax1[patch].axhline(y=sim.evolution[patch,2,-1],ls='--', label=r'$R_{inf}$')

         ax1[patch].legend()
         ax1[patch].set_xlabel("No. of iterations")
         ax1[patch].set_ylabel(r'$R_{inf}$')
         #ax1.set_title(r'$R(\infty)$ vs  $\mu/\beta$' )


    plt.show()

def interpretation_test():
    n = 2 #numero de parches

    #vector de parametros para una zona
    param = np.zeros(2)
    param[0] = beta = 0.67 # 0.67
    param[1] = gamma = 1./7. #1./7.

    #initial conditions para una zona
    y = np.zeros(3)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0


    p11=0.9
    p12=1.-p11
    p22=1.0
    p21=0.0
    P=np.array([[p11, p12],[p21,p22]])

    sirModel = SIR(n,params=[[beta,gamma],[beta+0.3,gamma]])
    sirModel2 = SIR(n,params=[[beta+0.3,gamma],[beta,gamma]])
    R0_1 = beta/gamma
    R0_2 = (beta+0.3)/gamma
    # Calculate Analitical Local final size
    R_infty = []
    R_infty_2 = []
    p_s = np.arange(0.,0.9,0.1)
    N = np.array([S,S])
    for p in p_s:
        P[0,0]=1.-p
        P[0,1]=p
        sirModel.set_network(MobilityNetwork(P))
        R_infty.append(sirModel.local_final_size(N,N))
        sirModel2.set_network(MobilityNetwork(P))
        R_infty_2.append(sirModel2.local_final_size(N,N))

    R_infty = np.array(R_infty)
    print(R_infty[:,0])
    plt.plot(p_s,R_infty[:,0], label="patch 1: R0=" + str(R0_1))
    plt.plot(p_s,R_infty[:,1], label="patch 2: R0=" + str(R0_2))
    plt.xlabel(r'$p_{1,2}$')
    plt.legend()

    plt.figure()
    R_infty_2 = np.array(R_infty_2)
    print(R_infty_2[:,0])
    plt.plot(p_s,R_infty_2[:,0], label="patch 1: R0=" + str(R0_2))
    plt.plot(p_s,R_infty_2[:,1], label="patch 2: R0=" + str(R0_1))
    plt.xlabel(r'$p_{1,2}$')
    plt.legend()
    plt.show()

    exit()

    sim = simulation(sirModel)
    sim.set_initial_conditions_all_patches(y)
    N = sim.node[:,0]
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(0,y) #Se establece esta condicion inicial en la zona 1
    sim.run()
    sim.plot_all()

def condition_test():
    n = 2 #numero de parches

    # parametros para una zona
    beta = 0.67 # 0.67
    gamma = 1./7. #0.1428

    #initial conditions para una zona
    N = np.array([35000, 35000]) #1500. 25000

    p11=0.9
    p12=1.-p11
    p22=1.0
    p21=0.0
    P=np.array([[p11, p12],[p21,p22]])
    sirModel = SIR(n,network=P)

    #simulacion
    y = np.zeros(3)
    S = y[0] = 35000 #1500. 25000
    I = y[1] = 0.0
    R = y[2] = 0.0
    sim = simulation(sirModel)
    sim.set_initial_conditions_all_patches(y)
    y[1]=y[1]+1. # Se agrega 1 infectado a las condiciones iniciales
    sim.set_initial_conditions_patch(0,y) #Se establece est

    Ro = []
    alphaAs = []
    R_inf = []
    R_inf_numerical = []
    betas = np.arange(0.02,0.8,0.1)
    for beta_ in betas:
        sirModel.set_patches_params(params=[[beta_,gamma],[beta,gamma]])
        cond, alphaA = sirModel.final_size_condition(N,N)
        R_inf.append(sirModel.local_final_size(N,N))
        Ro.append(beta_/gamma)
        alphaAs.append(alphaA)
        sim.run()
        R_inf_numerical.append(sim.evolution[:,2,-1])


    plt.plot(Ro,alphaAs, label=r'$|alpha| ||A||$')
    plt.axhline(y=1./np.exp(1.),ls='--', label=r'$1/e$')
    plt.xlabel(r'$R_o$')
    plt.legend()

    plt.figure()
    plt.plot(Ro,np.array(R_inf).T[0], label="Analitical R_inf_1")
    plt.plot(Ro,np.array(R_inf).T[1], label="Analitical R_inf_2")
    plt.plot(Ro,np.array(R_inf_numerical).T[0], label="Numerical R_inf_1")
    plt.plot(Ro,np.array(R_inf_numerical).T[1], label="Numerical R_inf_2")
    plt.xlabel(r'$R_o$')
    plt.legend()
    plt.show()


if (True):
    #index()
    #condition_test()
    interpretation_test()
