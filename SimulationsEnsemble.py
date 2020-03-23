from Simulation import simulation
import copy
from modelos import PLOSModel
from modelos import SIR
from MobilityNetwork import MobilityNetwork
import matplotlib.pyplot as plt
import matplotlib.lines as ml


class simulationsEnsemble:
	
	def __init__(self):
		self.simulations = []
		self.infected_average = None
		self.recoverd_average = None
		self.number_of_simulations = 0
		
	def add_simulation(self, simulation):
		self.simulations.append(copy.deepcopy(simulation))
		self.number_of_simulations += 1
	
	def run_simulation(self,i):
		self.simulations[i].run()
		
	def run_all_simulations(self):
		for i in range(self.number_of_simulations):
			self.run_simulation(i)
		self.average_infected()
		self.average_recovered()
			
	def average_infected(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.run_simulation(0)
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		infected_average = self.simulations[0].total_infected.copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				self.run_simulation(i)
			infected_average += self.simulations[i].total_infected
		self.infected_average = infected_average/(n)
		return infected_average.copy()
	
	def average_suseptible(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.run_simulation(0)
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		susceptible_average = self.simulations[0].total_susceptible.copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				self.run_simulation(i)
			susceptible_average += self.simulations[i].total_susceptible
		self.susceptible_average = susceptible_average/(n)
		return susceptible_average.copy()
		
	def average_recovered(self):
		n=self.number_of_simulations
		if (n>0 and self.simulations[0].runned_times == 0):
			self.run_simulation(0)
		elif (n==0):
			print("You have to add some simulations...")
			exit()
		recovered_average = self.simulations[0].total_recovered.copy()
		for i in range(1,n):
			if (self.simulations[i].runned_times == 0):
				self.run_simulation(i)
			recovered_average += self.simulations[i].total_recovered
		self.recovered_average = recovered_average/(n)
		return recovered_average.copy()
		
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
			self.run_simulation(n)
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

		
	def plot_infected_average(self, axes=None):
		if (self.infected_average is None ):
			self.average_infected()
		I_average = self.infected_average
		time = self.simulations[0].time
		if (axes is None):
			fig, ax = plt.subplots()
		else:
			ax = axes
		ax.set_xlabel('Time')
		ax.set_ylabel('Infected (simulation average)')
		ax.plot(time,I_average)
		#plt.show()
		return ax
		
	def plot_recovered_average(self, axes=None):
		if (self.recovered_average is None ):
			self.average_recovered()
		R_average = self.recovered_average
		time = self.simulations[0].time
		if (axes is None):
			fig, ax = plt.subplots()
		else:
			ax = axes
		ax.set_xlabel('Time')
		ax.set_ylabel('Recovered (simulation average)')
		ax.plot(time,R_average)
		#plt.show()
		return ax

if __name__ == '__main__':
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

#-----------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
	#numero de simulaciones en el ensemble
	number_of_simulations = 10
	
	#parameteres
	n = 6 #numero de parches
	p = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = [0]*2
	param[0] = beta1 = 1.5
	param[1] = gama1 = 1.
	#initial conditions para una zona
	x  = [0]*3
	S1 = x[0] = 2500.
	I1 = x[1] = 0.0 
	R1 = x[2] = 0.0

	
	sim = simulation() #se crea objeto simulacion
	sim.add_many_patches_parameters(n,param)  # se agregan n zonas con parametros dados en param
	sim.set_initial_conditions_all_patches(x) # se agregan las condiciones iniciales para todos las zonas
	x[1]=x[1]+1. # Se agrega 1 infecto a las condiciones iniciales
	sim.set_initial_conditions_patch(1,x) #Se establece esta condicion inicial en la zona 1
	sim.set_model(SIR) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	
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