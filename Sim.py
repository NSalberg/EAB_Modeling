import matplotlib.pyplot as plt
import numpy as np
import csv
import scipy.optimize
from scipy.integrate import odeint


weeks = range(0, 100)
deaths = []




# Parameter values
# 1/alpha = time spent as larvae 40-50 weeks
alpha = 1/45

# 1/delta = time spent as adult ~12 weeks
delta = 1/12

# number of eggs per adult
mu = 20

# time wasps spend as larvae
mu_prime = 1/45

# 1/delta_prime = time wasps spend as adults 4-6 weeks
delta_prime = 1/5

# We need to test these parameters
beta = 0.0000001
gamma = 35
K = 1000
params0 = np.array([alpha, delta, mu, mu_prime, delta_prime, beta, gamma, K])


y0 = [100, 0, 1, 0]

def sim(variables, t, params):
    B = variables[0]
    Bl = variables[1]
    Il = variables[2]
    W = variables[3]
    
    alpha = params[0]
    delta = params[1]
    mu = params[2]
    mu_prime = params[3]
    delta_prime = params[4]

    beta = params[5]
    gamma = params[6]
    k = params[7]


    dBdt = alpha * Bl - delta * B
    dBldt = mu * B - alpha * Bl - gamma * Bl**2 / (K**2 + Bl**2 + Il**2) - beta * Bl * W
    dIldt = beta * Bl * W - mu_prime * Il - gamma * Il**2 / (K**2 + Bl**2 + Il**2)
    dWdt = mu_prime * Il - delta_prime * W

    return [dBdt, dBldt, dIldt, dWdt]

t = np.linspace(weeks[0], weeks[-1], num=10000)

output = odeint(sim, y0, t, args=(params0,))







# Plot the results
f,EAB_and_EAB_larvae =  plt.subplots(1,1)
# add a title
EAB_and_EAB_larvae.set_title('EAB and Eab larvae vs time')


EAB_and_EAB_larvae.set_xlabel('Weeks')
EAB_and_EAB_larvae.set_ylabel('Population')

borer, = EAB_and_EAB_larvae.plot(t, output[:, 3], label='Borer')
borer_larvae, = EAB_and_EAB_larvae.plot(t, output[:, 1], label='Borer Larvae')

EAB_and_EAB_larvae.legend(handles=[borer_larvae, borer],)

#plot all the lines on a another graph
f,allLines =  plt.subplots(1,1)

# add a small title that describes the model params
allLines.set_title('EAB, Eab larvae, infected larvae, and wasps vs time')

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

plt.show()




