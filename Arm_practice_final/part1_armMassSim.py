import matplotlib.pyplot as plt
import armParam as P
from armAnimation import armAnimation
from dataPlotter import dataPlotter
from armDynamics import armDynamics
import numpy as np

# instantiate arm, controller, and reference classes
rod = armDynamics()

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = armAnimation()
# control gain calculations go here
th_eq = 0.0 *np.pi/180.0
tau_eq = P.m * P.g * P.ell / 2.0 * np.cos(th_eq)

t = P.t_start
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        u = tau_eq
        y = rod.update(u)
        t = t + P.Ts
    # update animation and data plots
    animation.update(rod.state)
    dataPlot.update(t, 0, rod.state, u)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
