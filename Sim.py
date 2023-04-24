import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.optimize
from scipy.integrate import odeint


weeks = range(0, 1000)
deaths = []




# Parameter values
# 1/alpha = time spent as larvae 40-50 weeks
alpha = 1/45

# 1/delta = time spent as adult ~12 weeks
delta = 1/12

# number of eggs per adult
mu = .5

# time wasps spend as larvae
mu_prime = 1/45

# 1/delta_prime = time wasps spend as adults 4-6 weeks
delta_prime = 1/5

# We need to test these parameters
beta = 0.0002
gamma = 35
K = 1000


params0 = np.array([alpha, delta, mu, mu_prime, delta_prime, beta, gamma, K])


y0 = [1000, 0, 0, 1000]

def sim(variables, t, params):
    # variables
    B = variables[0]
    Bl = variables[1]
    Il = variables[2]
    W = variables[3]


    # params 
    alpha = params[0]
    delta = params[1]
    mu = params[2]
    mu_prime = params[3]
    delta_prime = params[4]

    beta = params[5]
    gamma = params[6]
    K = params[7]


    dBdt = alpha * Bl - delta * B
    dBldt = mu * B - alpha * Bl - gamma * Bl**2 / (K**2 + Bl**2 + Il**2) - beta * Bl * W
    dIldt = beta * Bl * W - mu_prime * Il - gamma * Il**2 / (K**2 + Bl**2 + Il**2)
    dWdt = mu_prime * Il - delta_prime * W

    return [dBdt, dBldt, dIldt, dWdt]

t = np.linspace(weeks[0], weeks[-1], num=1000)

output = odeint(sim, y0, t, args=(params0,))


np.savetxt("output.txt",output, fmt="%.18f")



f, (EAB_and_Wasps, EAB_larvae_and_infected_larvae, allLines) = plt.subplots(3, 1, sharex=True)

# Plot the results
EAB_larvae_and_infected_larvae.set_title('Infected Larvae and EAB larvae vs time')
EAB_larvae_and_infected_larvae.set_xlabel('Weeks')
EAB_larvae_and_infected_larvae.set_ylabel('Population')

infected_larvae, = EAB_larvae_and_infected_larvae.plot(t, output[:, 2], label='Infected Larvae')
borer_larvae, = EAB_larvae_and_infected_larvae.plot(t, output[:, 1], label='Borer Larvae')

EAB_larvae_and_infected_larvae.legend(handles=[borer_larvae, infected_larvae],)


#plot all the lines on the same graph
# add a small title that describes the model params
allLines.set_title('EAB, EAB larvae, infected larvae, and wasps vs time')

"""
these are the annotations for the parameters

allLines.annotate('beta = ' + str(round(min[0],3) ) , xy=(.7, 150), xytext=(0, 0), textcoords='offset points') # type: ignore
allLines.annotate('gamma = ' + str(round(min[1],3)), xy=(.7, 125), xytext=(0, 0), textcoords='offset points') # type: ignore
allLines.annotate('delta = ' + str(round(min[2],3)), xy=(.7, 100), xytext=(0, 0), textcoords='offset points') # type: ignore
"""
allLines.set_xlabel('Weeks')
allLines.set_ylabel('Population')

B, = allLines.plot(t, output[:, 0], label='EAB')
Bl, = allLines.plot(t, output[:, 1], label='EAB larvae')
Il, = allLines.plot(t, output[:, 2], label='Infected larvae')
W, = allLines.plot(t, output[:, 3], label='Wasps')

allLines.legend(handles=[B, Bl, Il, W])


#plot EAB and wasps on a another graph
EAB_and_Wasps.set_title('EAB and Wasps vs time')

EAB_and_Wasps.set_xlabel('Weeks')
EAB_and_Wasps.set_ylabel('Population')

infected_larvae, = EAB_and_Wasps.plot(t, output[:, 0], label='EAB')
wasps, = EAB_and_Wasps.plot(t, output[:, 3], label='Wasps')

EAB_and_Wasps.legend(handles=[infected_larvae, wasps],)



plt.show()












