#  programs.py
#  
#  Copyright 2018 Jorge Velazquez Castro
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
import matplotlib.pyplot as plt
from Simulation import simulation
from modelos import PLOSModel
from MobilityNetwork import MobilityNetwork
from control_protocol import controlProtocol
from SimulationsEnsemble import simulationsEnsemble

def Homogeneus_Simple_Control():
	#parameteres

	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = np.zeros(5)
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.
	param[2] = alfa1 =5 
	param[3] = c1 =1000 
	param[4] = mu1 = 1./8.
	#initial conditions para una zona
	y = np.zeros(5)
	S1 = y[0] = 1500.
	I1 = y[1] = 0.0 
	R1 = y[2] = 0.0
	Sv1 = y[3] = 200 
	Iv1 = y[4] = 0.0
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	P.binomial(n,p,min_residential) #en este caso es una red binomial
	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_conectivity_network(P) # se agrega la red de conectividad P
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones inicials para todos las zonas
	y[1]=y[1]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	sim.set_simulation_time(80) #How many "days" to simulate 
	
	params=sim.parameters
	ControlSimple=controlProtocol(params,P)
	ControlSimple.set_observation_interval(7.)
	
	sim.set_control_protocol(ControlSimple)
	sim.run() #Se corre la simulacion
	sim.plot_all() # Se grafica I para la zona 0.
	sim.plot_total_infected()
	return sim
	
def Homogeneus_Without_Control():
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = np.zeros(5)
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.
	param[2] = alfa1 =5 
	param[3] = c1 =1000 
	param[4] = mu1 = 1./8.
	#initial conditions para una zona
	y = np.zeros(5)
	S1 = y[0] = 1500.
	I1 = y[1] = 0.0 
	R1 = y[2] = 0.0
	Sv1 = y[3] = 200 
	Iv1 = y[4] = 0.0
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	P.binomial(n,p,min_residential) #en este caso es una red binomial
	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_conectivity_network(P) # se agrega la red de conectividad P
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones inicials para todos las zonas
	y[1]=y[1]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	sim.set_simulation_time(80) #How many "days" to simulate 
	sim.run() #Se corre la simulacion
	sim.plot_all() # Se grafica I para la zona 0.
	sim.plot_total_infected()
	return sim
	
def Without_Control():
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = [0]*5
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.
	param[2] = alfa1 =5 
	param[3] = c1 =1000 
	param[4] = mu1 = 1./8.
	Capacities = [900.,500.,1500.,2000.,5000.,3000.]
	#initial conditions para una zona
	y = np.zeros(5)
	S1 = y[0] = 1500.
	I1 = y[1] = 0.0 
	R1 = y[2] = 0.0
	Sv1 = y[3] = 200 
	Iv1 = y[4] = 0.0
	
	sim = simulation() #se crea objeto simulacion
	for C in Capacities:
		param[3]=C
		sim.add_one_patch_parameters(list(param))  # se agregan n zonas con parametros dados en param
	
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones inicials para todos las zonas
	y[1]=y[1]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	P.binomial(n,p,min_residential) #en este caso es una red binomial
	sim.set_conectivity_network(P) # se agrega la red de conectividad P
	
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	sim.set_simulation_time(80) #How many "days" to simulate 
	sim.run() #Se corre la simulacion
	sim.plot_all() # Se grafica I para la zona 0.
	return sim
	
def simple_simulations_ensemble():
	#numero de simulaciones en el ensemble
	number_of_simulations = 8
	
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = [0]*5
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.
	param[2] = alfa1 =5 
	param[3] = c1 =1000 
	param[4] = mu1 = 1./8.
	#initial conditions para una zona
	y = [0]*5
	S1 = y[0] = 1500.
	I1 = y[1] = 0.0 
	R1 = y[2] = 0.0
	Sv1 = y[3] = 200 
	Iv1 = y[4] = 0.0
	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
	y[4]=y[4]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	
	ensemble = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones
	for i in range(number_of_simulations) :
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)
		ensemble.add_simulation(sim)

	ensemble.run_all_simulations() # run all simulations in the ensemble
	
	ensemble.plot_infected_average()
	
def MPI_simulations_ensemble():	#corre un ensamble de simulaciones en paralelo usando MPI. Correr con "mpiexec -n <number_of_process> python programs.py"
	#numero de simulaciones en el ensemble
	number_of_simulations = 8
	
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = [0]*5
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.
	param[2] = alfa1 =5 
	param[3] = c1 =1000 
	param[4] = mu1 = 1./8.
	#initial conditions para una zona
	y = [0]*5
	S1 = y[0] = 1500.
	I1 = y[1] = 0.0 
	R1 = y[2] = 0.0
	Sv1 = y[3] = 200 
	Iv1 = y[4] = 0.0
	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
	y[4]=y[4]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	
	ensemble = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones
	for i in range(number_of_simulations) :
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)
		ensemble.add_simulation(sim)

	ensemble.MPI_run_all() # run all simulations in the ensemble using MPI. Run using: "mpiexec -n <number_of_process> python programs.py"
	
	ensemble.plot_infected_average()
	
#def comparison_SimpleControl_NoControl():
	

if (__name__ == '__main__'):
	#Homogeneus_Simple_Control()
	#Homogeneus_Without_Control()
	#Without_Control()
	#simple_simulations_ensemble()
	#MPI_simulations_ensemble()  # requiere "mpiexec" 
