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
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from Simulation import simulation
from modelos import PLOSModel
from modelos import SIR
from MobilityNetwork import MobilityNetwork
from Control_protocol import controlProtocol
from SimulationsEnsemble import simulationsEnsemble
from Control_protocol import randomControl
from Control_protocol import finalEpidemicSizeControl
import copy

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
	#-----------------------------------------------------------------------------------------------------------------------

def SIR_Homogeneo():
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = np.zeros(2)
	param[0] = beta1 = 0.67
	param[1] = gama1 = 1./7.

	#initial conditions para una zona
	x = np.zeros(3)
	S1 = x[0] = 1500.
	I1 = x[1] = 0.0 
	R1 = x[2] = 0.0
	N  = [1500,1500,1500,1500,1500,1500]
	P=MobilityNetwork() #Se crea la red de mobilidad
	P.binomial(n,p,min_residential) #en este caso es una red binomial
	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_conectivity_network(P) # se agrega la red de conectividad P
	sim.set_initial_conditions_all_patches(x) # se agregan las condiciones inicials para todos las zonas
	x[1]=x[1]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,x) #Se establece esta condicion inicial en la zona 1
	sim.set_model(SIR) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	params=sim.parameters
	sim.set_simulation_time(80) #How many "days" to simulate 
	FinalSizeControl=finalEpidemicSizeControl(params,P,N)
	print(FinalSizeControl.R_inf)
	print(FinalSizeControl.Thetas)
	sim.set_control_protocol(FinalSizeControl)
	sim.run() #Se corre la simulacion
	sim.plot_all() # Se grafica I para la zona 0.
	sim.plot_total_infected()
	#sim.control_protocol.Itotal(N)
	return sim


	#-----------------------------------------------------------------------------------------------------------------------
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
	plt.show()
	
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
	plt.show()
	
def comparison_SimpleControl_NoControl():
	#numero de simulaciones en el ensemble
	number_of_simulations = 100
	
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
	
	#Se genera el ensamble de simulaciones sin control
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
	y[4]=y[4]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	
	P=MobilityNetwork() #Se crea la red de mobilidad
	
	##Sin Control
	ensembleNoControl = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones (contendra simulaciones sin control)
	for i in range(number_of_simulations) :
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)
		ensembleNoControl.add_simulation(sim)

	ensembleNoControl.run_all_simulations() # run all simulations in the ensemble
	
	##Control Simple cada 7 dias:
	#Se adjunta el control y se genera el ensamble de simulaciones con control
	ensembleSimpleControl = simulationsEnsemble() #Se crea objeto ensable de simulaciones
	for i in range(number_of_simulations) :
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)  #se establece la red de conectividad generada para ser usada en la simulacion
		params = sim.parameters
		ControlSimple=controlProtocol(params,P)  #Se crea objeto control (necesita los parametros del modelo)
		ControlSimple.set_observation_interval(7.) #Se establecen los dias para ajustar el control
		sim.set_control_protocol(ControlSimple) #Se agrega el control especificado a la simulacion
		ensembleSimpleControl.add_simulation(sim) #Se agrega la simulacion al conjunto de simulaciones a ejecutar (al ensamble)
		
	ensembleSimpleControl.run_all_simulations() 
	
	##Control Simple cada dia:
	##ensembleSimpleControlDaily = copy.deepcopy(ensembleSimpleControl)
	##for sim1 in ensembleSimpleControlDaily.simulations :
	##	sim1.control_protocol.set_observation_interval(1.)
	
	ensembleSimpleControlDaily = simulationsEnsemble() #Se crea objeto ensable de simulaciones
	for i in range(number_of_simulations) :
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)  #se establece la red de conectividad generada para ser usada en la simulacion
		params = sim.parameters
		ControlSimple=controlProtocol(params,P)  #Se crea objeto control (necesita los parametros del modelo)
		ControlSimple.set_observation_interval(1.) #Se establecen los dias para ajustar el control
		sim.set_control_protocol(ControlSimple) #Se agrega el control especificado a la simulacion
		ensembleSimpleControlDaily.add_simulation(sim) #Se agrega la simulacion al conjunto de simulaciones a ejecutar (al ensamble)
	
	ensembleSimpleControlDaily.run_all_simulations()
	
	
	####Plot total infected of both cases ####
	
	axes = ensembleNoControl.plot_infected_average()
	axes.lines[0].set_label('No Control')
	axes1 = ensembleSimpleControl.plot_infected_average(axes)
	axes1.lines[1].set_label('Simple Control(7 dias)')
	axes2 = ensembleSimpleControlDaily.plot_infected_average(axes1)
	axes2.lines[2].set_label('Simple Control(diario)')
	
	plt.legend()
	plt.show()
	
def comparison_SimpleControl_NoControl_SameNetwork():
	#numero de simulaciones en el ensemble
	number_of_simulations = 100
	
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # elementos de la diagonal de la matriz de mobilidad mayor a este numero
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
	
	ensembleNoControl = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones (contendra simulaciones sin control)
	ensembleSimpleControl = simulationsEnsemble()
	ensembleSimpleControlDaily = simulationsEnsemble()
	for i in range(number_of_simulations) :
		#Se genera la simulacion base
		sim = simulation() #se crea objeto simulacion
		sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
		y[1]=0
		sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
		y[1]=1. # Se agrega 1 infectado a las condiciones iniciales
		sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
		sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
		P=MobilityNetwork() #Se crea la red de mobilidad
		P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
		sim.set_conectivity_network(P)  #se establece la red de conectividad generada para ser usada en la simulacion
		#No control
		ensembleNoControl.add_simulation(sim) ##No control
		#Control Semanal
		params = sim.parameters
		ControlSimple=controlProtocol(params,P)  #Se crea objeto control (necesita los parametros del modelo)
		ControlSimple.set_observation_interval(7.) #Se establecen los dias para ajustar el control
		sim.set_control_protocol(ControlSimple) #Se agrega el control especificado a la simulacion
		ensembleSimpleControl.add_simulation(sim) #Se agrega la simulacion al conjunto de simulaciones a ejecutar (al ensamble)
		#Control diario:
		ControlSimple=controlProtocol(params,P) 
		ControlSimple.set_observation_interval(1.)
		sim.set_control_protocol(ControlSimple)
		ensembleSimpleControlDaily.add_simulation(sim)
		
	ensembleNoControl.run_all_simulations()
	ensembleSimpleControl.run_all_simulations()
	ensembleSimpleControlDaily.run_all_simulations()
	
	####Plot total infected of all cases ####
	
	axes = ensembleNoControl.plot_infected_average()
	axes.lines[0].set_label('No Control')
	axes1 = ensembleSimpleControl.plot_infected_average(axes)
	axes1.lines[1].set_label('Simple Control(7 dias)')
	axes2 = ensembleSimpleControlDaily.plot_infected_average(axes1)
	axes2.lines[2].set_label('Simple Control(diario)')
	
	plt.legend()
	plt.show()
	
def MPI_control_comparison():
	#numero de simulaciones en el ensemble
	number_of_simulations = 240 #8
	
	#parameteres
	n = 30 #numero de parches
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
	
	import mpi4py
	from mpi4py import MPI
	
	comm = MPI.COMM_WORLD
	size = comm.Get_size()
	rank = comm.Get_rank()
	
	
	if rank == 0 :
		ensembleNoControl = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones (contendra simulaciones sin control)
	if rank == 1 :
		ensembleSimpleControl = simulationsEnsemble()
	if rank == 2:
		ensembleSimpleControlDaily = simulationsEnsemble()
	if rank == 3:
		ensembleRandomControl = simulationsEnsemble()
	
	for i in range(number_of_simulations) :
		#Se genera la simulacion base
		sim = simulation() #se crea objeto simulacion
		sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
		y[1]=0
		sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
		y[1]=1. # Se agrega 1 infectado a las condiciones iniciales
		sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
		sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
		sim.set_simulation_time(200)  #Se establece el tiempo de simulacion.
		if rank == 0:
			P=MobilityNetwork() #Se crea la red de mobilidad
			#P.binomial(n,p,min_residential) #para cada simulacion se genera una nueva red binomial
			#P.barabsi_albert(n,2,min_residential) # argumentos (n=No_de_nodos,m=No_de_enlaces que se incertan al incertar un nodo,min_residential)
			P.random(n,5*n,min_residential) # argumentos (n=No_de_nodos,m=No_de_enlaces,min_residential), red muy conectada
			#P.random(n,n,min_residential) # red aleatoria poco conectada
		else:
			P=None
		P=comm.bcast(P,root=0)
		sim.set_conectivity_network(P)  #se establece la red de conectividad generada para ser usada en la simulacion
		if rank == 0:
			#No control
			ensembleNoControl.add_simulation(sim) ##No control
		if rank == 1:
			#Control Semanal
			params = sim.parameters
			ControlSimple=controlProtocol(params,P)  #Se crea objeto control (necesita los parametros del modelo)
			ControlSimple.set_observation_interval(7.) #Se establecen los dias para ajustar el control
			sim.set_control_protocol(ControlSimple) #Se agrega el control especificado a la simulacion
			ensembleSimpleControl.add_simulation(sim) #Se agrega la simulacion al conjunto de simulaciones a ejecutar (al ensamble)
		if rank == 2:
			#Control diario:
			params = sim.parameters
			ControlSimple=controlProtocol(params,P) 
			ControlSimple.set_observation_interval(1.)
			sim.set_control_protocol(ControlSimple)
			ensembleSimpleControlDaily.add_simulation(sim)
		if rank == 3:
			#Random control diario:
			params = sim.parameters
			ControlRandom=randomControl(params,P) 
			ControlRandom.set_observation_interval(7.)
			sim.set_control_protocol(ControlRandom)
			ensembleRandomControl.add_simulation(sim)



	if rank==0:
		print("Antes de hacer las sim...para " + str(rank))
		ensembleNoControl.run_all_simulations()
		print("Simulaciones hechas para " + str(rank))
	if rank==1:
		print("Antes de hacer las sim...para " + str(rank))
		ensembleSimpleControl.run_all_simulations()
		print("Simulaciones hechas para " + str(rank))
	if rank==2:
		print("Antes de hacer las sim...para " + str(rank))
		ensembleSimpleControlDaily.run_all_simulations()
		print("Simulaciones hechas para " + str(rank))
	if rank==3:
		print("Antes de hacer las sim...para " + str(rank))
		ensembleRandomControl.run_all_simulations()
		print("Simulaciones hechas para " + str(rank))

	comm.barrier()
	####Plot total infected of all cases ####
	if rank == 1:
		print("enviando simulaciones de " + str(rank))
		comm.send(ensembleSimpleControl, dest=0)
	elif rank == 0:
		ensembleSimpleControl = comm.recv(source=1)
		
	if rank == 2:
		print("enviando simulaciones de " + str(rank))
		comm.send(ensembleSimpleControlDaily, dest=0)
	elif rank == 0:
		ensembleSimpleControlDaily = comm.recv(source=2)
		
	if rank == 3:
		print("enviando simulaciones de " + str(rank))
		comm.send(ensembleRandomControl, dest=0)
	elif rank == 0:
		ensembleRandomControl = comm.recv(source=3)
	
	if(rank == 0):
		print("Simulaciones recibidas en " + str(rank))
	else:
		print("simulacion ha sido enviada desde " + str(rank))
	
	comm.barrier()
	if rank == 0:
		axes = ensembleNoControl.plot_infected_average()
		axes.lines[0].set_label('No Control')
		axes1 = ensembleSimpleControl.plot_infected_average(axes)
		axes1.lines[1].set_label('Simple Control(7 dias)')
		axes2 = ensembleSimpleControlDaily.plot_infected_average(axes1)
		axes2.lines[2].set_label('Simple Control(diario)')
		axes3 = ensembleRandomControl.plot_infected_average(axes2)
		axes3.lines[3].set_label('Ramdom Control (diario)')
		
		plt.legend()
		plt.show()
		
		axesR = ensembleNoControl.plot_recovered_average()
		axesR.lines[0].set_label('No Control')
		axesR1 = ensembleSimpleControl.plot_recovered_average(axesR)
		axesR1.lines[1].set_label('Simple Control(7 dias)')
		axesR2 = ensembleSimpleControlDaily.plot_recovered_average(axesR1)
		axesR2.lines[2].set_label('Simple Control(diario)')
		axesR3 = ensembleRandomControl.plot_recovered_average(axesR2)
		axesR3.lines[3].set_label('Ramdom Control (diario)')
		
		plt.legend()
		plt.show()
	
		

if (__name__ == '__main__'):
	#Homogeneus_Simple_Control()
	#Homogeneus_Without_Control()
	#Without_Control()
	#simple_simulations_ensemble()
	#MPI_simulations_ensemble()  # requiere MPI instalado : "mpiexec" Run using: "mpiexec -n <number_of_process> python programs.py" (n=8)
	SIR_Homogeneo()
	#comparison_SimpleControl_NoControl()
	#comparison_SimpleControl_NoControl_SameNetwork()
	#MPI_control_comparison() # Run using: "mpiexec -n 4 python programs.py"
