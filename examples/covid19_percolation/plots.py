import numpy as np
import matplotlib.pyplot as plt



max_infected = np.load("max_infected_car_nc.npy")
total_infected = np.load("total_infected_car_nc.npy")

p_nodes = np.load('p_nodes.npy')
p_edges = np.load('p_edges.npy')

np.savetxt("max_infected_car_nc.csv", max_infected, delimiter=",")
np.savetxt("total_infected_car_nc.csv", total_infected, delimiter=",")
np.savetxt("p_nodes.csv", p_nodes, delimiter=",")
np.savetxt("p_edges.csv", p_edges, delimiter=",")


x, y = np.meshgrid(p_nodes, p_edges)
z_test= x*y

fig, ax = plt.subplots()
z = max_infected[:-1,:-1]
z=z.T
z_min = z.min()
z_max = z.max()

print(z_test.shape)
print(z.shape)


c = ax.pcolormesh(x, y, z, cmap='hot', vmin=z_min, vmax=z_max)
ax.set_title('Peak of Infection, Mun+Car')
#ax.set_title('Size of infection, Mun+Car')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_xlabel(r'$q_n$')
ax.set_ylabel(r'$q_a$')
fig.colorbar(c, ax=ax)

#######
fig, ax = plt.subplots()

z = total_infected[:-1,:-1]
z=z.T
z_min = z.min()
z_max = z.max()

print(z_test.shape)
print(z.shape)


c = ax.pcolormesh(x, y, z, cmap='hot', vmin=z_min, vmax=z_max)
ax.set_title('Size of infection, Mun+Car')
# set the limits of the plot to the limits of the data
ax.axis([x.min(), x.max(), y.min(), y.max()])
ax.set_xlabel(r'$q_n$')
ax.set_ylabel(r'$q_a$')
fig.colorbar(c, ax=ax)



plt.show()

exit()
plt.imshow(max_infected)
plt.show()
