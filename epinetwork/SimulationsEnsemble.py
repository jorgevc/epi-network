from .Simulation import simulation
import copy
from .MobilityNetwork import MobilityNetwork
import matplotlib.pyplot as plt
import matplotlib.lines as ml
import numpy as np


class simulationsEnsemble:
	__slots__=('simulations','susceptible_average','infected_average',\
	'recovered_average','number_of_simulations','max_infected','evolution_average')

	def __init__(self):
		self.simulations = []
		self.susceptible_average = None
		self.infected_average = None
		self.recovered_average = None
		self.number_of_simulations = 0
		self.max_infected = None
		self.evolution_average = None

	def add_simulation(self, simulation):
		self.simulations.append(copy.deepcopy(simulation))
		if simulation.time == None:
			self.simulations[-1].time = np.arange(0.,simulation.simulation_time,1.)
		self.number_of_simulations += 1

	def run_all_simulations(self):
		for i in range(self.number_of_simulations):
			self.simulations[i].run()
		self.average_infected()
		self.average_recovered()

	def get_evolution_average(self,mode="relative"):
		if mode=="relative":
			evols_list = [x.get_evolution_totals()/np.sum(x.evolution_totals[:,0]) for x in self.simulations]
		else:
			evols_list = [x.get_evolution_totals() for x in self.simulations]
		evols = np.array(evols_list)
		self.evolution_average = np.sum(evols,axis=0)/len(evols_list)
		return self.evolution_average, self.simulations[0].time

	def average_infected(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.simulations[0].run()
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		infected_var = self.simulations[0].infected_variable
		infected_average = self.simulations[0].evolution_totals[infected_var,:].copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				print("Warning: Running in non multiprocessing ")
				self.simulations[i].run()
			infected_average += self.simulations[i].evolution_totals[infected_var,:]
		self.infected_average = infected_average/(n)
		return infected_average.copy()

	def average_suseptible(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.run_simulation(0)
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		susceptible = self.simulations[0].susceptible_variable
		susceptible_average = self.simulations[0].evolution_totals[susceptible,:].copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				self.run_simulation(i)
			susceptible_average += self.simulations[i].evolution_totals[susceptible,:]
		self.susceptible_average = susceptible_average/(n)
		return susceptible_average.copy()

	def average_recovered(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.run_simulation(0)
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		recovered=self.simulations[0].recovered_variable
		recovered_average = self.simulations[0].evolution_totals[recovered,:].copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				self.run_simulation(i)
			recovered_average += self.simulations[i].evolution_totals[recovered,:]
		self.recovered_average = recovered_average/(n)
		return recovered_average.copy()

	def calculate_max_infected(self):
		if (self.infected_average is None):
			self.average_infected()
		self.max_infected=max(self.infected_average)
		return self.max_infected.copy()

	def MPI_run_all(self):
		import  mpi4py
		from mpi4py import MPI

		comm = MPI.COMM_WORLD
		size = comm.Get_size()
		rank = comm.Get_rank()

		number_process_simulations = [int(self.number_of_simulations/size)]*size
		number_process_simulations[0] += self.number_of_simulations % size
		process_offset = [0]*size
		process_simulations = []

		for i in range(1,size):
			process_offset[i]=process_offset[i-1] + number_process_simulations[i-1]

		for i in range(number_process_simulations[rank]):
			n=process_offset[rank] + i
			self.simulations[n].run()
			process_simulations.append(self.simulations[n])

		comm.barrier()
		recvbuf = comm.gather(process_simulations, root=0)
		if rank==0:
			self.simulations = [item for sublist in recvbuf for item in sublist]
			self.average_infected()
			self.average_recovered()
			MPI.Finalize()
		else:
			MPI.Finalize()
			exit()

	def multiprocessing_run_all(self, No_of_process=8):
		from multiprocessing import Pool
		s=self.number_of_simulations
		with Pool(processes=No_of_process) as pool:
			doing = [pool.apply_async(self.simulations[i].run, ()) for i in range(s)]
			simulations = [res.get() for res in doing]
		self.simulations = simulations

	def save_average(self,file_name):
		try:
			import cPickle as pickle
		except ModuleNotFoundError:
			import pickle

		ensemble = simulationsEnsemble()
		#sim = simulation()
		#for i in range(self.number_of_simulations):
		#	sim = copy.deepcopy(self.simulations[i])
		#	sim.evolution=[]
		#	ensemble.add_simulation(sim)

		ensemble.infected_average = self.infected_average
		ensemble.recovered_average = self.recovered_average

		with open(file_name + ".SimulationsEnsemble", "wb") as f:
			pickle.dump(ensemble, f, pickle.HIGHEST_PROTOCOL)

	def plot_infected_average(self, axes=None, relative=False):
		if (self.infected_average is None ):
			self.average_infected()
		I_average = self.infected_average
		time = self.simulations[0].time
		if (axes is None):
			fig, ax = plt.subplots()
		else:
			ax = axes
		ax.set_xlabel('Time')
		if (relative==True):
			if (self.recovered_average is None ):
				self.average_recovered()
			R_ini=self.recovered_average[0]
			if (self.susceptible_average is None):
				self.average_suseptible()
			S_ini=self.susceptible_average[0]
			I_ini=self.infected_average[0]
			ax.set_ylabel('% of Infected Population')
			ax.plot(time,(I_average/(R_ini + S_ini+I_ini))*100.)
		else:
			ax.set_ylabel('Infected (simulation average)')
			ax.plot(time,I_average)
		#plt.show()
		return ax

	def plot_recovered_average(self, axes=None, relative=False):
		if (self.recovered_average is None ):
			self.average_recovered()
		R_average = self.recovered_average
		time = self.simulations[0].time
		if (axes is None):
			fig, ax = plt.subplots()
		else:
			ax = axes
		ax.set_xlabel('Time')
		if (relative==True):
			if (self.infected_average is None ):
				self.average_infected()
			I_ini=self.infected_average[0]
			if (self.susceptible_average is None):
				self.average_suseptible()
			S_ini=self.susceptible_average[0]
			R_ini=self.recovered_average[0]
			ax.set_ylabel('Recovered (% population)')
			ax.plot(time,(R_average/(R_ini + S_ini+I_ini))*100.)
		else:
			ax.set_ylabel('Recovered (simulation average)')
			ax.plot(time,R_average)
		#plt.show()
		return ax

if __name__ == '__main__':
	from .modelos import PLOSModel
	#numero de simulaciones en el ensemble
	number_of_simulations = 10

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
	ensemble.plot_recovered_average()
	plt.show()
