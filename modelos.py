#  modelos.py
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
#  se espera recibir (yv:state, t:time, param:array[], p:matrix[Network],n:No_patches)
import numpy as np

def PLOSModel(yv,t,param,p,n,control):
# regresa la derivada del sistema de ecuaciones diferenciales
   y = yv.reshape(n,5)

   Sh = np.zeros(n)
   Ih = np.zeros(n)
   Rh = np.zeros(n)
   Sv = np.zeros(n)
   Iv = np.zeros(n)
   Nh = np.zeros(n)
   Nv = np.zeros(n)
   W = np.zeros(n)
   F = np.zeros(n)

   beta = np.zeros(n)
   gamma = np.zeros(n)
   alfa = np.zeros(n)
   c = np.zeros(n)
   mu = np.zeros(n)

   j=0
   z=0


   for i in range (0,n):
      Sh[i] = y[i][0]
      Ih[i] = y[i][1]
      Rh[i] = y[i][2]
      Sv[i] = y[i][3]
      Iv[i] = y[i][4]

      Nh[i] = Sh[i] + Ih[i] + Rh[i]
      Nv[i] = Sv[i] + Iv[i]


   if(control.observation_time(t)):
         control.observe(y,t)

   control.calculate_control(t)

   for i in range (0,n):
      control_params = control.get_control(i)
      beta[i] = param[i][0]*(1. - control_params[0])
      gamma[i] = param[i][1] # + control_params[1]
      alfa[i] = param[i][2]
      c[i] = param[i][3]
      mu[i] = param[i][4]

      W[i] = sum(Nh*p[:,i])
      F[i] = Iv[i] / W[i]

   for i in range (0,n):

      dS = -beta[i]*Sh[i]*sum(p[i,:]*F)
      dI = beta[i]*Sh[i]*sum(p[i,:]*F)-gamma[i]*Ih[i]
      dR = gamma[i]*Ih[i]
      dSv = alfa[i]*Nv[i]*(1-Nv[i]/c[i])-beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Sv[i]
      dIv = beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Iv[i]

      if i==0:
          b = np.array([dS,dI,dR,dSv,dIv])
      else:
          a = np.array([dS,dI,dR,dSv,dIv])
          b = np.concatenate((b,a),axis=0)

   return b

def SEAIR(yv,t,param,p,n,control):
    y = yv.reshape(n,5)

    S = np.zeros(n)
    E = np.zeros(n)
    A = np.zeros(n)
    I = np.zeros(n)
    R = np.zeros(n)
    N = np.zeros(n)
    W = np.zeros(n)
    F = np.zeros(n)
    f = np.zeros(n)

    #params:
    beta = np.zeros(n)
    sigma = np.zeros(n)
    gamma = np.zeros(n)
    mu = np.zeros(n)
    m = np.zeros(n)

    S = y[:,0]
    E = y[:,1]
    A = y[:,2]
    I = y[:,3]
    R = y[:,4]

    N = S + E + A + I + R
    #params
    if (t>15.):
        tau=55.
        Ctr = np.exp(-(t-15.)/tau)
    else:
        Ctr = 1.

    beta = np.array(param)[:,0]*Ctr
    sigma = np.array(param)[:,1]
    gamma = np.array(param)[:,2]
    m = np.array(param)[:,3]
    eta= np.array(param)[:,4]

    F_I=I.dot(p)
    F_A=A.dot(p)
    F_E=E.dot(p)
    W=N.dot(p)

    dS = - beta * S * p.dot( (F_E + F_A + eta*F_I)/W  )
    dE = beta * S * p.dot( (F_E + F_A + eta*F_I)/W  ) - sigma*E
    dA = m*sigma*E - gamma*A
    dI = (1.-m)*sigma*E-gamma*I
    dR = gamma*(A + I)

    return np.array([dS, dE, dA, dI, dR]).T.flatten()

def NPatch2(yv,param,p,n):
   y = yv.reshape(n,5)
# regresa la solucion del sistema de ecuaciones diferenciales con la condicion (beta[i]*Sh[i]*sum(p[i,:]*F)/sigma[i])>1 para cada i

   Sh = np.zeros(n)
   Ih = np.zeros(n)
   Rh = np.zeros(n)
   Sv = np.zeros(n)
   Iv = np.zeros(n)
   Nh = np.zeros(n)
   Nv = np.zeros(n)
   W = np.zeros(n)
   F = np.zeros(n)

   beta = np.zeros(n)
   sigma = np.zeros(n)
   alfa = np.zeros(n)
   c = np.zeros(n)
   mu = np.zeros(n)

   j=0
   z=0

   for i in range (0,n):
      Sh[i] = y[i][0]
      Ih[i] = y[i][1]
      Rh[i] = y[i][2]
      Sv[i] = y[i][3]
      Iv[i] = y[i][4]

      Nh[i] = Sh[i] + Ih[i] + Rh[i]
      Nv[i] = Sv[i] + Iv[i]

      beta[i] = param[i][0]
      sigma[i] = param[i][1]
      alfa[i] = param[i][2]
      c[i] = param[i][3]
      mu[i] = param[i][4]

   for i in range (0,n):
      W[i] = sum(Nh*p[:,i])
      F[i] = Iv[i] / W[i]


   for i in range (0,n):

      if  (beta[i]*Sh[i]*sum(p[i,:]*F)/sigma[i])>1:

         dS = -beta[i]*Sh[i]*sum(p[i,:]*F)
         dI = beta[i]*Sh[i]*sum(p[i,:]*F)-sigma[i]*Ih[i]
         dR = sigma[i]*Ih[i]
         dSv = alfa[i]*Nv[i]*(1-Nv[i]/c[i])-beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Sv[i]
         dIv = beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Iv[i]
      else:

         dS = -beta[i]*Sh[i]*sum(p[i,:]*F)
         dI = -sigma[i]*Ih[i]
         dR = sigma[i]*Ih[i]
         dSv = alfa[i]*Nv[i]*(1-Nv[i]/c[i])-beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Sv[i]
         dIv = beta[i]*(Sv[i]/W[i])*sum(Ih*p[:,i])-mu[i]*Iv[i]


      if i==0:
          b = np.array([dS,dI,dR,dSv,dIv])
      else:
          a = np.array([dS,dI,dR,dSv,dIv])
          b = np.concatenate((b,a),axis=0)
   return b
