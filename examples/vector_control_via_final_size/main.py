



def main():
    #parameteres

	n = 6 #numero de parches
	b = 0.5 # parametro de la red binomial
	min_residential = 0.8 # diagonal de la matriz de mobilidad mayor a este numero
	#vector de parametros para una zona
	param = np.zeros(5)
	param[0] = beta_h = 0.67
	param[1] = gamma = 1./7.
	param[2] = beta_v =5
	param[3] = mu_v = 1./8.
	#initial conditions para una zona
	y = np.zeros(5)
	S = y[0] = 1500.
	I = y[1] = 0.0
	R = y[2] = 0.0
	V = y[3] = 200
	W = y[4] = 0.0


    P = MobilityNetwork()
    P.binomial(n,b,min_residential)

    vectorModel = VectorBorne(network=P)
    vectorModel.set_patches_params(param)

    sim = simulation(vectorModel)
    sim.set_initial_conditions_all_patches(y)
    y[1]=y[1]+1. # Se agrega 1 infectado mosquito a las condiciones iniciales
	sim.set_initial_conditions_patch(1,y) #Se establece esta condicion inicial en la zona 1
	sim.set_simulation_time(80) #How many "days" to simulate
    sim.run() #Se corre la simulacion

	sim.plot_total_infected()
	sim.plot_total_recovered()




if __name__ == '__main__':
    main()
