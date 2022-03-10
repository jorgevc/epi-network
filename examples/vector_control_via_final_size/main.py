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
from epinetwork.Control_protocol import randomControl, IndexBasedControl

import matplotlib
import matplotlib.pyplot as plt

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

    vectorModel = VectorBorne(params=param, network=P, control=IndexBasedControl)

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



if (True):
    main()
