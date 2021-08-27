import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline
import csaps

peak_covid=np.load('infection_peak_covid_n.npy')
peak_normal=np.load('infection_peak_normal_n.npy')
size_covid=np.load('infection_size_covid_n.npy')
size_normal=np.load('infection_size_normal_n.npy')
links=np.load('n_links_n.npy')

n=30 #100.
S=1000.
Norm=100./(n*S)

x=links
# Peak
# y=(1. - peak_covid/peak_normal)*100.
# spl= UnivariateSpline(x,y, s=1490000)
# #cs = CubicSpline(links,peak_normal-peak_covid)
#
# plt.plot(links/n,peak_normal*Norm,'o',label="base line mobility")
# plt.plot(links/n,peak_covid*Norm,'+',label="reduced mobility")
#
# plt.xlabel("links/node")
# plt.ylabel("% of infected population")
# plt.legend()
#
# plt.figure()
#
#
# xs = np.linspace(15,90,50)
#
#
# plt.plot(links/n, (1. - peak_covid/peak_normal)*100.,'o',label="simulations data")
# plt.plot(xs/n,spl(xs),label="spline interpolation")
# plt.xlabel("links per node")
# plt.ylabel( "% of reduced infections")
# plt.legend()
#
# plt.show()

#size
y=(1. - size_covid/size_normal)*100.
spl= UnivariateSpline(x,y, s=1490000)
#cs = CubicSpline(links,peak_normal-peak_covid)

plt.plot(links/n,size_normal*Norm,'o',label="base line mobility")
plt.plot(links/n,size_covid*Norm,'+',label="reduced mobility")

plt.xlabel("links/node")
plt.ylabel("% of infected population")
plt.legend()

plt.figure()


xs = np.linspace(15,90,50)


plt.plot(links/n, (1. - size_covid/size_normal)*100.,'o',label="simulations data")
plt.plot(xs/n,spl(xs),label="spline interpolation")
plt.xlabel("links per node")
plt.ylabel( "% of reduced infections")
plt.legend()

plt.show()
