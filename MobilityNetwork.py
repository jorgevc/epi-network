#  MobilityNetwork.py
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
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random as rand


class MobilityNetwork:
	
	def __init__(self):
		self.matrix = None
		self.network = None
		
	def weight_edges(self, Net, min_residential):
		"""
		Give weights to the edges of the of a networkx network, it also
		adds the corresponding weight to the diagonal in order to 
		represent a mobility network where the rows are normalized. 
		"""
		
		n=nx.number_of_nodes(Net)
		for i in range (0,n):
			Net.add_edge(i,i)
		
		for u,v,d in Net.edges(data=True):
			if u == v:
				d['weight']=rand.uniform(min_residential,1)
			else:
				d['weight']=rand.uniform(0,min_residential)
		
		B=nx.to_numpy_matrix(Net,weight='weight')
		
		for i in range (0,n):
			suma=0
			for j in range (0,n):
				suma=suma+B[i,j]
			B[i,:]=B[i,:]/suma
			
		NewNet=nx.from_numpy_matrix(B)

		self.matrix = np.squeeze(np.asarray(B))
		self.network = NewNet
		
		return B
		

	def binomial(self,n,p,min_residential): 
		"""
		Se crea un matriz en el cual se muestra la interaccion de los parches
		
		Parameters
		--------
		n : int
			Numero de nodos en la red
		p : float
			Parametro de la matriz binomial
		min_residential: float
			Minima fraccion de tiempo que los individuos estan en sus propios parches.
			Es decir la diagonal tendra numeros mayores o iguales a este valor
		
		Returns
		------
		numpy matriz of size : n x n
		"""
		
		
		G=nx.binomial_graph(n,p)
		self.weight_edges(G, min_residential)

	def binomial2(self,n,p2,min_residential): 
		"""
		Se crea un matriz en el cual se muestra la interaccion de los parches
		
		Parameters
		--------
		n : int
			Numero de nodos en la red
		p : float
			Parametro de la matriz binomial
		min_residential: float
			Minima fraccion de tiempo que los individuos estan en sus propios parches.
			Es decir la diagonal tendra numeros mayores o iguales a este valor
		
		Returns
		------
		numpy matriz of size : n x n
		"""
		
		
		G=nx.binomial_graph(n,p2)
		self.weight_edges(G, min_residential)

		
	def barabsi_albert(self,n,m,min_residential):
		"""
		Se crea un matriz en el cual se muestra la interaccion de los 
		parches con el algoritmo de Barabasi albert (red de mundo pequeno)
		
		Parameters
		--------
		n : int
			Numero de nodos en la red
		m : float
			Number of edges to attach from a new node to exiting nodes
		min_residential: float
			Minima fraccion de tiempo que los individuos estan en sus propios parches.
			Es decir la diagonal tendra numeros mayores o iguales a este valor
		
		Returns
		------
		numpy matriz of size : n x n
		"""
		
		G=nx.barabasi_albert_graph(n,m)
		self.weight_edges(G, min_residential)

	def random(self,n,m,min_residential):
		
		G=nx.gnm_random_graph(n,m)
		self.weight_edges(G, min_residential)

		
		
	def draw_network():
		nx.draw_networkx(self.network)
		plt.show()
		
