import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.interpolate import CubicSpline
import csaps

peak_covid=np.load('infection_peak_covid.npy')
peak_normal=np.load('infection_peak_normal.npy')
size_covid=np.load('infection_size_covid.npy')
size_normal=np.load('infection_size_normal.npy')
links=np.load('n_links.npy')

n=100.
S=1000.
Norm=100./(n*S)

x=links
y=peak_normal-peak_covid
spl= UnivariateSpline(x,y, s=1490000)
#cs = CubicSpline(links,peak_normal-peak_covid)

plt.plot(links/n,peak_normal*Norm,'o',label="base line mobility")
plt.plot(links/n,peak_covid*Norm,'+',label="stay at home mobility")

plt.xlabel("links/node")
plt.ylabel("% of infected population")
plt.legend()

plt.figure()


xs = np.linspace(30,400,100)


plt.plot(links/n,Norm*(peak_normal-peak_covid),'o',label="simulations data")
plt.plot(xs/n,Norm*spl(xs),label="spline interpolation")
plt.xlabel("links per node")
plt.ylabel( "% of reduced infections")
plt.legend()

plt.show()
