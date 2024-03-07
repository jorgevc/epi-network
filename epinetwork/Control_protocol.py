#  control_protocol.py
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
import random as rnd

class Protocol:
	def __init__(self,*,update_interval=1.,model=None):
		self.observation = None
		self.last_observation_time = None
		self.observation_interval = update_interval
		self.control = None
		self.last_control_time = None
		self.control_interval = update_interval
		self.model = model

	def get_control(self,y,t):
		if(self.observation_time(t)):
			self.observation=y.copy()
			self.last_observation_time = t
		if(self.last_control_time is None):
			self.calculate_control()
			self.last_control_time = t
		elif( self.control_interval is not None):
			if(t-self.last_control_time > self.control_interval):
				self.calculate_control()
				self.last_control_time = t
		return self.control

	def calculate_control(self):
		pass

	def observe(self,y,t):
		self.observation=y.copy()
		self.last_observation_time=t

	def observation_time(self,t):
		if (self.last_observation_time is None):
			return True
		elif ((self.observation_interval is not None)):
			if ((t - self.last_observation_time) > self.observation_interval):
				return True
		else:
			return True

	def set_observation_interval(self, dt):
		self.observation_interval = dt
		if(self.control_interval==None):
			self.control_interval = dt

	def set_control_interval(self, dt):
		self.control_interval = dt
		if(self.observation_interval == None):
			self.observation_interval=dt

	def set_model(self, model):
		self.model = model


class controlProtocol(Protocol):
	def __init__(self,params,P_network):
		self.number_of_patches = len(params)
		self.observations = [[0]*5]*self.number_of_patches
		self.params = params
		self.P_network = P_network
		self.last_observation_time = None
		self.observation_interval = 1.
		self.control = [[0.,0.]]*self.number_of_patches
		self.TRh=[] #human Transmission index
		self.VR=[]  #Vulnerability index (human)
		self.TRv=[] #vector Transmission idex
		self.last_control_time = 0.
		self.control_interval = 1.

	def observe(self,y,t):
		self.observations=y.copy()
		self.last_observation_time=t
		self.calculate_indices(t)

	def recalculate_control(self):
		TRh = self.TRh[-1][1]
		TRh_max = max(TRh)
		max_transmition_patch = TRh.index(TRh_max)
		#if( TRh_max > 2.):
		self.control = [[0.,0.]] * self.number_of_patches
		self.control[max_transmition_patch] = [0.5,0.5]

	def calculate_control(self,t):
		if(t-self.last_control_time > self.control_interval):
			self.recalculate_control()
		return self.control

	def calculate_indices(self,t):
		n=self.number_of_patches
		Sh=[]
		Ih=[]
		Rh=[]
		Sv=[]
		Iv=[]
		beta=[]
		gamma=[]
		mu=[]
		W=[]
		Nh=[]

		for i in range(n): 		# i corre sobre los parches
			Sh.append(self.observations[i][0])
			Ih.append(self.observations[i][1])
			Rh.append(self.observations[i][2])
			Sv.append(self.observations[i][3])
			Iv.append(self.observations[i][4])
			beta.append(self.params[i][0])
			gamma.append(self.params[i][1])
			mu.append(self.params[i][4])
			Nh.append(Sh[i] + Ih[i] + Rh[i])

		Rv=np.zeros((n,n))
		Rh=np.zeros((n,n))
		p=self.P_network.matrix
		for k in range(n):
			W.append(np.dot(Nh,p[:,k]))
			for j in range(n):
				Rv[k,j]=(beta[j]*Sh[j]*p[j,k]*Iv[k])/(mu[k]*W[k]) #Secondary human infections of residentes from patch j produced by infected vectors in k (See ec. 20 of ref [1])
				Rh[j,k]=beta[k]*Ih[j]*p[j,k]*Sv[k]/(gamma[j]*W[k]) #Secondary vector cases in patch k caused by infected residents of patch i

		self.TRh.append((t,[(sum(np.dot(Rh,Rv)[i,:])) for i in range(n)])) #secondary infections generated by a single individual of patch i
		self.VR.append((t,[(sum(np.dot(Rh,Rv)[:,j])) for j in range(n)])) #vulnerability of patch j
		self.TRv.append((t, [sum(Rv[i,:]) for i in range(n)])) #vector transmission index (verificar)
		return self.TRh[-1][1]

	def plot_TRindex(self, i):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('TR_i')
		plt.plot(self.time,self.evolution[i,1,:])
		plt.show()


class noControl(Protocol):
	def __init__(self,model=None):
		super().__init__(update_interval=None,model=model)

	def calculate_control(self):
		self.control = np.zeros((self.model.number_of_patches,1))


class RandomControl(Protocol):
	def calculate_control(self):
		random_patch = rnd.randint(0,self.model.number_of_patches-1)
		self.control = np.zeros((self.model.number_of_patches,1))
		self.control[random_patch] = np.array([1.])

class IndexBasedControl(Protocol):
	def __init__(self,*,strategy='only_max',cost=None,update_interval=1.,model=None):
		super().__init__(update_interval=update_interval,model=model)
		self.cost = cost
		self.cost_to_ctrl_K = 1. #When cost per node is 1 the ctrl is 0.5
		if (self.model != None):
			self.num_contrl_patches = model.number_of_patches
		if strategy.lower() == 'only_max' :
			self.calculate_control = self.calculate_only_max
		elif strategy.lower() == 'equal_distributed' :
			self.calculate_control = self.calculate_equal_distributed
		elif strategy.lower() == 'cost_prop' :
			self.calculate_control = self.calculate_cost_prop_distributed
		elif strategy.lower() == 'control_prop' :
			self.calculate_control = self.calculate_prop_control


	def set_model(self, model):
		self.model = model
		self.num_contrl_patches = model.number_of_patches

	def set_cost(self,cost):
		self.cost = cost

	def set_number_of_controlled_patches(self,num):
		self.num_contrl_patches = num

	def cost_to_ctrl(self,cost):
		K=self.cost_to_ctrl_K
		return cost/(K+cost)

	def calculate_equal_distributed(self):
		self.indices = self.model.get_indices(self.observation)
		
		n=self.num_contrl_patches
		if n<len(self.indices):
			ctrl_patches = np.argpartition(self.indices,-n)[-n:]
		else:
			ctrl_patches = np.arange(len(self.indices))

		self.control = np.zeros((self.model.number_of_patches,1))
		self.control[ctrl_patches,0]=self.cost_to_ctrl(self.cost/n)

	def calculate_cost_prop_distributed(self):
		self.indices = self.model.get_indices(self.observation)
		n=self.num_contrl_patches
		if n<len(self.indices):
			ctrl_patches = np.argpartition(self.indices,-n)[-n:]
		else:
			ctrl_patches = np.arange(len(self.indices))

		cost_distribution = (self.indices[ctrl_patches]*self.cost)/np.sum(self.indices[ctrl_patches])
		self.control = np.zeros((self.model.number_of_patches,1))
		self.control[ctrl_patches,0]=self.cost_to_ctrl(cost_distribution)

	def calculate_prop_control(self):
		epochs = 100
		learning_rate = 0.01
		max_error = 0.01
		self.indices = self.model.get_indices(self.observation)
		n=self.num_contrl_patches
		if n<len(self.indices):
			ctrl_patches = np.argpartition(self.indices,-n)[-n:]
		else:
			ctrl_patches = np.arange(len(self.indices))

		def C(sigma, ctrl_patches):
			S=self.indices[ctrl_patches]*sigma #TODO Avoid S gd 1
			K=self.cost_to_ctrl_K
			return K*np.sum(S/(1.-S))
		
		it=0
		sigma=0.
		max_index = np.max(self.indices[ctrl_patches])
		learning_rate = learning_rate/max_index
		err = C(sigma,ctrl_patches) - self.cost
		while (abs(err)>max_error and it < epochs):
			sigma=sigma - learning_rate*err
			sigma = 0.95/max_index if sigma*max_index > 0.95 else sigma
			err = C(sigma,ctrl_patches) - self.cost
			it+=1

		print("Cost error = ", C(sigma,ctrl_patches)-self.cost)
		
		self.control = np.zeros((self.model.number_of_patches,1))
		self.control[ctrl_patches,0]=sigma*self.indices[ctrl_patches]


	def calculate_only_max(self):
		self.indices = self.model.get_indices(self.observation)
		patch_index = np.argmax(self.indices)
		self.control = np.zeros((self.model.number_of_patches,1))
		self.control[patch_index,0] = self.cost_to_ctrl(self.cost)




