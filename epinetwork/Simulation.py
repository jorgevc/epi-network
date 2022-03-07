#  Simulation.py
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
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt
from .MobilityNetwork import MobilityNetwork
from .Control_protocol import controlProtocol
from .Control_protocol import noControl
import copy
import inspect

class simulation:
	"""Class simulation represents a single simulation of a epidemic

    Attributes
    ----------
    The following attributes need to be set in order to run the simulation:
        No_patches : int
          Number of patches (zones) in the network.
        parameters : list[list[float]]
          List of the parameter's list for each patch
        model : function object
          function representing the ODE of the dinamical system
        P : 2dim numpy array
          Conectivity matrix of the network
        node : list of numpy arrays
          List of initial conditions of the patches(nodes).
        simulation_time : int
          final time in the corresponding units of the model parameters of the simulation
        control_protocol : controlProtocol object
          object representing the control protocol to be implemented (Default non control)

    The following attributes are calculated and set after the simulation is runned with the 'run' procedure:
        time : numpy.array
          vector of times where state of the system was computed
        evolution : 3 dimentional numpy array
          numpy array[node][equation][step time] giving the estate of the variable (equation) at a time step
        total_infected : numpy.array[i])
          Total infected of the simulation at time steep i (sum of infected of all patches) 2nd ecuation is assumed to be infected
        total_recovered : numpy.array[i]
          Total recovered of the simulation at time steep i (sum of recovered of all patches) 3rd equation is assumed to be recovered
        runned_times : int
          Number of times the simulation has been run.

    Methods
    -------
    add_one_patch_parameters(parameters):
       Adds one patch (sector) to the simulation with list of parameters given as the argument
    """

	__slots__=('No_patches','parameters','node','time','model','P','evolution',\
	'control_protocol','simulation_time','total_infected','total_susceptible',\
	'total_recovered','runned_times','patch_dimention','infected_variable',\
	'recovered_variable')
	def __init__(self):
		#parameteres
		self.No_patches = 0 #numero de parches
		self.parameters = []
		self.node = []
		self.time = []
		self.model = None
		self.P = None
		self.evolution = None
		self.control_protocol = noControl()
		self.simulation_time = 100
		self.total_infected = None
		self.total_susceptible = None
		self.total_recovered = None
		self.runned_times = 0
		self.patch_dimention = 0 #dimention of de system of a single patch
		self.infected_variable = 1
		self.recovered_variable = 2

	def reset(self):
		self.__init__()

	def add_one_patch_parameters(self,parameters):
		"""Adds one patch (sector) to the simulation with list of parameters given as the argument

		Parameters
		----------
		parameters : list or 1dim numpy array
		   list of patch parameters

		Returns
		-------
		none

		"""

		self.parameters.append(list(parameters))
		self.node.append([])
		self.No_patches +=1

	def add_many_patches_parameters(self,No_patches, parameters):
		"""Adds many similar patches to the simulation

		Parameters
		----------
		No_parches : int
		   number of patches to be added to the simulation
		parameters : list of 1dim numpy array
		   list of parameters of the patch. All patches will have this same parameters

		Returns
		-------
		none

		"""

		for i in range(No_patches):
			self.add_one_patch_parameters(parameters)

	def set_initial_conditions_patch(self,i,Initial):
		"""Sets the initial conditions of a specific patch

		Parameters
		----------
		i : int
			Number of patch
		Initial : list or numpy array
			list of values for the initial conditions of a single patch

		Returns
		-------
		none

		"""

		self.node[i] = np.array(Initial)
		if self.patch_dimention == 0:
			self.patch_dimention = self.node[i].shape[0]

	def set_initial_conditions_all_patches(self,Initial):
		"""Sets the same initial conditions for all the patches

		Parameters
		----------
		Initial : list or 1 dim numpy array
		   List of initial conditions for a patch. This same conditions will be set for all patches

		Returns
		-------
		none

		"""

		for i in range(self.No_patches):
			self.set_initial_conditions_patch(i,Initial)

	def set_conectivity_network(self,P):
		self.P=copy.deepcopy(P)

	def set_model(self,model):
		if isinstance(model,Model):
			self.model = model
		else:
			n=self.No_patches
			p = self.P.matrix
			control = self.control_protocol
			parameters = self.parameters
			def system(t,estate):
				return model(estate,t,parameters,p,n,control)
			self.model.system=system

	def set_simulation_time(self,time):
		self.simulation_time=time

	def run(self):
		dim = self.patch_dimention
		self.time = np.linspace(0.0,self.simulation_time,self.simulation_time*10)

		#system(t,yv):
		system = self.model.system

		initial = np.array(self.node).flatten()
		#solution = odeint(system,initial,self.time)
		#self.evolution = [[solution[:,node*dim+column] for column in range(dim)] for node in range(n)]
		solution = np.array(solve_ivp(system, [min(self.time), max(self.time)], initial, t_eval=self.time).y)
		self.evolution = solution.reshape(n,dim,-1)

		self.calculate_total_infected()
		self.calculate_total_susceptible()
		self.calculate_total_recovered()
		self.runned_times += 1
		return self

	def set_control_protocol(self,control_protocol):
		self.control_protocol = copy.deepcopy(control_protocol)

	def graficar(self,i):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('Infectados')
		plt.plot(self.time,self.evolution[i,self.infected_variable,:])
		plt.show()

	def plot_all(self):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('Infectados')
		for i in range(self.No_patches) :
			plt.plot(self.time,self.evolution[i,self.infected_variable,:])
		plt.show()

	def calculate_total_infected(self):
		infected = self.infected_variable
		self.total_infected = self.evolution[0,infected,:]  #asumed 2nd equation is infected
		for i in range(1,self.No_patches):
			self.total_infected += self.evolution[i,infected,:]

	def calculate_total_susceptible(self):
		self.total_susceptible = self.evolution[0,0,:]  #asumed 1st equation is susceptible
		for i in range(1,self.No_patches):
			self.total_susceptible += self.evolution[i,0,:]

	def calculate_total_recovered(self):
		recovered = self.recovered_variable
		self.total_recovered = self.evolution[0,recovered,:]  #asumed 3er ecuation is recovered
		for i in range(1,self.No_patches):
			self.total_recovered += self.evolution[i,recovered,:]

	def plot_total_infected(self):
		if (self.total_infected is None ):
			self.calculate_total_infected()
		plt.figure()
		plt.xlabel('Time')
		plt.ylabel('Total Infected')
		line = copy.copy(plt.plot(self.time,self.total_infected))
		plt.show()
		return line

	def plot_total_recovered(self):
		if (self.total_recovered is None ):
			self.calculate_total_recovered()
		plt.figure()
		plt.xlabel('Time')
		plt.ylabel('Total Recovered')
		line = copy.copy(plt.plot(self.time,self.total_recovered))
		plt.show()
		return line

if __name__ == '__main__':
	from .modelos import PLOSModel
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
	sim.set_initial_conditions_all_patches(y) # se agregan las condiciones iniciales para todos las zonas
	y[4]=y[4]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion

	params=sim.parameters
	ControlSimple=controlProtocol(params,P)
	ControlSimple.set_observation_interval(5.)

	sim.set_control_protocol(ControlSimple)
	sim.run() #Se corre la simulacion

	sim.plot_total_infected()
	sim.plot_total_recovered()

	#sim.plot_all() # Se grafica I para la zona 0.
