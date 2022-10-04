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
from .MobilityNetwork import MobilityNetwork
from .Control_protocol import noControl, Protocol

class Model:
    def __init__(self,n=0, *, params = None , network=MobilityNetwork, control = noControl()):
     #Number of patches, parameter of a patch, mobility network, conctrol class
        self.number_of_variables = 0
        self.number_of_patches = n
        self.set_network(network)
        self.set_protocol(control)
        self.params = np.array(params)

    def system(self,t,yv):
        pass

    def set_protocol(self,control):
        if isinstance(control,Protocol):
            self.control = control
            self.control.model = self
        else:
            self.control = control(model = self)

    def set_network(self, network):
        if isinstance(network,type):
            self.p = network(self.number_of_patches)
        else:
            self.p = network
            if ( network.network.number_of_nodes() > self.number_of_patches ):
                self.number_of_patches = network.network.number_of_nodes()
                print("Caution : Number of nodes has been changed due to size of network")

class VectorBorne(Model):
    def __init__(self,n=1,*,params = np.full((4),None) , network=MobilityNetwork, control = noControl()):
         #Number of patches, parameter of a patch, mobility network, conctrol class
        super().__init__(n=n,params=params,network=network,control=control)
        self.number_of_variables = 5
        self.final_size_presition = 1.
        self.final_size_max_iterations = 200
        self.beta_h = np.full((self.number_of_patches),params[0])
        self.gamma = np.full((self.number_of_patches),params[1])
        self.beta_v = np.full((self.number_of_patches),params[2])
        self.mu_v = np.full((self.number_of_patches),params[3])

    def system(self,t,yv):
        y = yv.reshape(self.number_of_patches,self.number_of_variables)

        S = y[:,0]
        I = y[:,1]
        R = y[:,2]
        V = y[:,3]
        W = y[:,4]

        N = S + I + R
        Nv = V + W

        P = N.dot(self.p.matrix)
        F_I = self.p.matrix.dot(W/P)
        F_V = I.dot(self.p.matrix)

        u = self.control.get_control(y,t)

        dS = - self.beta_h * S * F_I
        dI =  self.beta_h * S * F_I - self.gamma * I
        dR = self.gamma * I
        dV = self.mu_v * Nv - self.beta_v * V * F_V/P - self.mu_v * (1.+ u[:,0]) * V
        dW = self.beta_v * V * F_V/P - self.mu_v * (1.+ u[:,0]) * W

        return np.array([dS, dI, dR, dV, dW ]).T.flatten()

    def  set_patches_params(self,params,No_patches=None):
        if (No_patches == None) :
            n = self.number_of_patches
        else :
            n = No_patches
            if (self.number_of_patches != n):
                self.number_of_patches = n
                print("Caution : Number of patches has changed")

        self.beta_h = np.full((n),params[0])
        self.gamma = np.full((n),params[1])
        self.beta_v = np.full((n),params[2])
        self.mu_v = np.full((n),params[3])


    def local_final_size(self,N,S,Nv):
        P = N.dot(self.p.matrix)
        b = self.p.matrix.dot(1./P)
        a = (self.beta_h * self.beta_v * Nv * b)/(self.gamma * self.mu_v * P)
        def f(x):
            tau = a * x.dot(self.p.matrix)
            return N - S*np.exp(-tau)

        err=1000
        it=0
        R_infty = f(np.ones(len(N)))
        while (err>self.final_size_presition or it < self.final_size_max_iterations):
            R_infty_next = f(R_infty)
            err = np.max(np.abs(R_infty_next - R_infty))
            it += 1
            R_infty = R_infty_next

        return R_infty

    def get_indices(self, y):

        S = y[:,0]
        I = y[:,1]
        R = y[:,2]
        V = y[:,3]
        W = y[:,4]

        N = S + I + R
        Nv = V + W

        return self.local_final_size(N,N,Nv)

class SIR(Model):
    def __init__(self,n=1, *, params = np.full((2),None) , network=MobilityNetwork, control = noControl()):
     #Number of patches, parameter of a patch, mobility network, conctrol class
        super().__init__(n=n,params=params,network=network,control=control)
        self.number_of_variables = 3
        self.final_size_presition = 0.01
        self.final_size_max_iterations = 200
        self.beta = None
        self.gamma = None
        self.set_patches_params()

    def system(self,t,yv):
        y = yv.reshape(self.number_of_patches,self.number_of_variables)

        S = y[:,0]
        I = y[:,1]
        R = y[:,2]
        N = S + I + R

        F_I=I.dot(self.p.matrix)
        W=N.dot(self.p.matrix)

        dS = -S*self.p.matrix.dot(self.beta*F_I/W)
        dI = S*self.p.matrix.dot(self.beta*F_I/W) - self.gamma*I
        dR = self.gamma*I

        return np.array([dS, dI, dR]).T.flatten()

    def set_patches_params(self,params=None,No_patches=None):
        if (No_patches == None) :
            n = self.number_of_patches
        else :
            n = No_patches
            if (self.number_of_patches != n):
                self.number_of_patches = n
                print("Caution : Number of patches has changed")
        if (params!=None):
            self.params=np.array(params)

        if (self.params.ndim == 1 and len(self.params)==2):
            self.beta = np.full((n),self.params[0])
            self.gamma = np.full((n),self.params[1])
        elif (self.params.ndim == 2 and len(self.params)==n):
            self.beta = self.params.T[0]
            self.gamma = self.params.T[1]
        else:
            print("Error: params have not been given or params for each patch is not the same of number of patches")


    def local_final_size(self,N,S):
        P = N.dot(self.p.matrix)
        b = np.multiply(self.p.matrix,self.beta/P)
        def f(x):
            tau = b.dot( (x/self.gamma).dot(self.p.matrix) )
            return N - S*np.exp(-tau)

        err=1000
        it=0
        R_infty = f(np.ones(len(N)))
        while (err>self.final_size_presition and it < self.final_size_max_iterations):
            R_infty_next = f(R_infty)
            err = np.max(np.abs(R_infty_next - R_infty))
            it += 1
            R_infty = R_infty_next

        return R_infty

    def get_indices(self, y):

        S = y[:,0]
        I = y[:,1]
        R = y[:,2]

        N = S + I + R

        return self.local_final_size(N,N)

    def final_size_condition(self, N, S_0):
        P = N.dot(self.p.matrix)
        b = np.multiply(self.p.matrix,self.beta/P).dot(self.p.matrix.T) # b_ij = sum_k (b_k*p_ik*p_jk)/P_k
        A = np.multiply(b,1./self.gamma)

        alpha = S_0*np.exp(-b.dot(N/self.gamma))

        alpha_A = np.norm(alpha)*np.norm(A)
        if (alpha_A < 1/np.exp(1)):
            condition = True
        else:
            condition = False

        return satisfy, alpha_A
