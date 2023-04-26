import matplotlib.pyplot as plt
import slopedMassParam as P
from slopedMassAnimation import slopedMassAnimation
from dataPlotter import dataPlotter
from slopedMassDynamics import slopedMassDynamics
import numpy as np

# instantiate arm, controller, and reference classes
slopedMass = slopedMassDynamics()

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = slopedMassAnimation()
# control gain calculations go here
z_eq = 0.0 
F_eq = P.k1 * z_eq + P.k2 * z_eq**3 - P.m * P.g * np.sin(np.pi/4)

t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        u = F_eq
        y = slopedMass.update(u)
        t = t + P.Ts
    # update animation and data plots
    animation.update(slopedMass.state)
    dataPlot.update(t, 0, slopedMass.state, u)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
