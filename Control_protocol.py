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

class controlProtocol:
	def __init__(self,params,P_network):
		self.number_of_patches = len(params)
		self.observations = [[0]*5]*self.number_of_patches
		self.params = params
		self.P_network = P_network
		self.last_observation_time = None
		self.observation_interval = 1.
		self.control = [[0.,0.]]*self.number_of_patches
		self.TRh=[]
		self.last_control_time = 0.
		self.control_interval = 1.

	def observe(self,y,t):
		self.observations=y
		self.last_observation_time=t
		self.calculate_indices(t)
		
	def observation_time(self,t):
		if(self.last_observation_time == None or (t - self.last_observation_time) > self.observation_interval):
			return True
		else:
			return False
		
	def recalculate_control(self):
		TRh = self.TRh[-1][1]
		TRh_max = max(TRh)
		max_transmition_patch = TRh.index(TRh_max)
		if( TRh_max > 2.):
			self.control = [[0.,0.]] * self.number_of_patches
			self.control[max_transmition_patch] = [0.2,0.2]
			
	def calculate_control(self,t):
		if(t-self.last_control_time > self.control_interval):
			self.recalculate_control()
		
	def get_control(self,i):
		return self.control[i]
		
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
		for i in range(n):
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
				Rv[k,j]=(beta[j]*Sh[j]*p[j,k]*Iv[k])/(mu[k]*W[k])
				Rh[j,k]=beta[k]*Ih[j]*p[j,k]*Sv[k]/(gamma[j]*W[k])
		
		self.TRh.append((t,[(sum(np.dot(Rh,Rv)[i,:])) for i in range(n)]))
		return self.TRh[-1][1]

	def set_observation_interval(self, dt):
		self.observation_interval = dt
		if(self.control_interval==None):
			self.control_interval = dt
		
	def set_control_interval(self, dt):
		self.control_interval = dt
		if(self.observation_interval == None):
			self.observation_interval=dt
		
	def plot_TRindex(self, i):
		plt.figure()
		plt.xlabel('Tiempo')
		plt.ylabel('TR_i')
		plt.plot(self.time,self.evolution[i][1][:])
		plt.show()

class noControl:
	def __init__(self,default = 0):
		self.default_control= [default,default]
		
	def observe(self,y,t):
		pass

	def observation_time(self,t):
		return False
		
	def calculate_control(self,t):
		pass
		
	def get_control(self,i):
		return self.default_control
