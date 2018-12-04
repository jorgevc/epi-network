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
import numpy as np
import matplotlib.pyplot as plt
from modelos import PLOSModel
from MobilityNetwork import MobilityNetwork
from control_protocol import controlProtocol
from control_protocol import noControl
import copy

class simulation:
	
	def __init__(self):
		#parameteres
		self.No_patches = 0 #numero de parches
		self.p = 0.
		self.min_residential_fraction = 1.
		self.parameters = []
		self.node = []
		self.time = []
		self.model = None
		self.P = None
		self.evolution = []
		self.control_protocol = noControl()
		self.simulation_time = 100
		
	def add_one_patch_parameters(self,parameters):
		self.parameters.append(list(parameters))
		self.node.append([])
		self.No_patches +=1 
		
	def add_many_patches_parameters(self,No_patches, parameters):
		for i in range(No_patches):
			self.add_one_patch_parameters(parameters)
			
	def set_initial_conditions_patch(self,i,Initial):
		self.node[i] = list(Initial)
		
	def set_initial_conditions_all_patches(self,Initial):
		for i in range(self.No_patches):
			self.set_initial_conditions_patch(i,Initial)
			
	def set_conectivity_network(self,P):
		self.P=copy.deepcopy(P)
		
	def set_model(self,model):
		self.model=model
	
	def set_simulation_time(self,time):
		self.simulation_time=time

	def run(self):
		n=self.No_patches
		p = self.P.matrix
		control = self.control_protocol
		parameters = self.parameters
		self.time = np.linspace(0.0,self.simulation_time,self.simulation_time*10)

		def system(estate,t):
			return self.model(estate,t,parameters,p,n,control)

		initial = np.array(self.node).flatten()
		solution = odeint(system,initial,self.time)
		self.evolution = [[solution[:,node*5+column] for column in range(5)] for node in range(n)]
		
	def set_control_protocol(self,control_protocol):
		self.control_protocol = copy.deepcopy(control_protocol)

	def graficar(self,i):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('Infectados')
		plt.plot(self.time,self.evolution[i][1][:])
		plt.show()
		
	def plot_all(self):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('Infectados')
		for i in range(self.No_patches) :
			plt.plot(self.time,self.evolution[i][1][:])
		plt.show()

if __name__ == '__main__':
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
	y[4]=y[4]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_model(PLOSModel) #Se establece el modelo usado en PLOS como modelo para hacer la simulacion
	
	params=sim.parameters
	ControlSimple=controlProtocol(params,P)
	ControlSimple.set_observation_interval(5.)
	
	sim.set_control_protocol(ControlSimple)
	sim.run() #Se corre la simulacion
	sim.plot_all() # Se grafica I para la zona 0.
