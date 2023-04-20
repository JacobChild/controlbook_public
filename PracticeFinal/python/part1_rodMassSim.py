import matplotlib.pyplot as plt
import rodMassParam as P
from rodMassAnimation import rodMassAnimation
from dataPlotter import dataPlotter
from rodMassDynamics import rodMassDynamics
import numpy as np

# instantiate arm, controller, and reference classes
rodMass = rodMassDynamics()

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = rodMassAnimation()
# control gain calculations go here
th_eq = 0.0 * np.pi/180.0
tau_eq =  P.m*P.g*P.ell * np.cos(th_eq) + P.k1 * th_eq + P.k2 * th_eq**3

t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        u = tau_eq
        y = rodMass.update(u)
        t = t + P.Ts
    # update animation and data plots
    animation.update(rodMass.state)
    dataPlot.update(t, th_eq, rodMass.state, u)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
