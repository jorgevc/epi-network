import numpy as np
import matplotlib
#matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt

from epinetwork.Simulation import simulation
from epinetwork.modelos import SEAIRD
from epinetwork.MobilityNetwork import MobilityNetwork
from epinetwork.Control_protocol import controlProtocol
from epinetwork.SimulationsEnsemble import simulationsEnsemble
from epinetwork.Control_protocol import randomControl
import copy
import random
import time

def Initial_infected_patch(N):
    N_total=sum(N)
    rand=random.uniform(0.,N_total)
    S=0.0
    for i in range(len(N)):
        S+=N[i]
        if rand <= S :
            return i

def delete_random_edges(graph,k=1):
    '''
     delete one edge of its current edge if del_orig is True.
    :param graph: networkx graph
    :param k: number of edges to delete
    :return: networkx graph
    '''
    edges = list(graph.edges)
    # random edges
    chosen_edges = random.sample(edges, k)

    for chosen_edge in chosen_edges:
        graph.remove_edge(chosen_edge[0], chosen_edge[1])

    return graph

def remove_random_nodes(G,k=1):
    nodes = list(G.nodes())
    chosen_nodes = random.sample(nodes,k)
    for chosen_node in chosen_nodes:
        G.remove_node(chosen_node)

    return G


if __name__ == '__main__':
    P = MobilityNetwork()
    B = P.from_adjacency_matrix("Datos/ma2.dat",0.8)
    n = B.shape[0] #number of patches
    no_edges = len(list(P.network.edges()))
    Population = np.loadtxt("Datos/poblacion.dat")

    P.add_nodes_property("Population",Population)

    #number_of_simulations = 5*226
    number_of_simulations = 500 #226

    units=50. # las delta de paso son realmente son 1/20
    p_nodes_max = int(0.6*units)
    p_edges_max = int(0.8*units)

    max_infected = np.zeros((p_nodes_max+1,p_edges_max+1))
    total_infected = np.zeros((p_nodes_max+1,p_edges_max+1))
    np.save('p_nodes.npy',np.arange(p_nodes_max)/units)
    np.save('p_edges.npy',np.arange(p_edges_max)/units)
    #parameteres
    min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
    #vector de parametros para una zona
    Ro=(2. + 4.5)/2.

    param = [0]*6

    param[1] = sigma = (0.2 + 0.33)/2.
    param[2] = gamma = (0.055 + 0.2)/2.
    param[3] = mu = 0.02
    param[4] = m = (0.2 + 0.8)/2.
    param[5] = eta = 0.2
    param[0] = beta = (Ro*gamma)/(eta*(1.-m)+m+gamma/sigma)
    #initial conditions para una zona
    y = [0]*6
    S = y[0] = 0.
    E = y[1] = 0.0
    A = y[2] = 0.0
    I = y[3] = 0.0
    R = y[4] = 0.0
    D= y[5] = 0.0
    elapsed_time=0.
    evaluations=0
    for p_nodes in range(p_nodes_max+1):
        for p_edges in range(p_edges_max+1):
            for i in range(4):
                start_time= time.time()
                n_to_delete= int(n*p_nodes/units)
                e_to_delete= int(no_edges*(p_edges/units))

                G=copy.deepcopy(P.network)
                G=delete_random_edges(G,e_to_delete)
                G=remove_random_nodes(G,n_to_delete)

                Network=MobilityNetwork()
                Network.from_netowrkx(G)
                no_patches=Network.matrix.shape[0]

                sim = simulation() #se crea objeto simulacion
                sim.add_many_patches_parameters(no_patches,param)  # se agregan n zonas con parametros dados en param
                sim.set_model(SEAIRD) #Se establece el modelo para hacer la simulacion
                sim.infected_variable = 3
                sim.recovered_variable = 4
                sim.set_simulation_time(360)
                New_population = []
                for idx, node in enumerate(Network.network.nodes(data="Population")):
                    y[0]=node[1]
                    sim.set_initial_conditions_patch(idx,y) #Todos suseptibles condicion inicial
                    New_population.append(node[1])

                ensemble = simulationsEnsemble() # se crea el objeto ensamble_de_simulaciones
                point1_time=time.time()
                for i in range(number_of_simulations) :
                    infected_patch=Initial_infected_patch(np.array(New_population))
                    y[1]=1. # Se agrega 1 expuesto a las condiciones iniciales
                    sim.set_initial_conditions_patch(infected_patch,y) #Se establece esta condicion inicial en la zona
                    Network.weight_edges(min_residential) #para cada simulacion se genera nueva mobilidad
                    sim.set_conectivity_network(Network)
                    ensemble.add_simulation(sim)
                    y[1]=0.0
                    sim.set_initial_conditions_patch(infected_patch,y) #se elimina el infectado para la siguiente simulacion
                point2_time=time.time()
                print("presimulation time ", point2_time - point1_time)
                ensemble.multiprocessing_run_all(12) #12 is the number of processes
                ensemble.calculate_max_infected()
                max_infected[p_nodes,p_edges] += (ensemble.max_infected/4.)
                ensemble.average_suseptible()
                R = ensemble.susceptible_average[0] - ensemble.susceptible_average[-1]
                total_infected[p_nodes,p_edges] += (R/4.)

                evaluations+=1.
                evals_to_go=p_nodes_max*p_edges_max*4 - evaluations
                elapsed_time+=time.time() - start_time
                average_time=elapsed_time/evaluations
                eta = average_time*evals_to_go
                print("-----", evaluations ," evaluations done-----")
                print( evals_to_go," evaluatons to go")
                print("aprox", eta/60. , "minutes left")
                print("or ", eta/3600, "hours.")

        if(p_nodes % 10 == 0):
            np.save("max_infected_car_nc.npy",max_infected)
            np.save("total_infected_car_nc.npy",total_infected)

    np.save("max_infected_car_nc.npy",max_infected)
    np.save("total_infected_car_nc.npy",total_infected)
