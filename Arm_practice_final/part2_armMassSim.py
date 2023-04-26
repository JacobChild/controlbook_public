import matplotlib.pyplot as plt
import numpy as np
import armParam as P
from signalGenerator import signalGenerator
from armAnimation import armAnimation
from dataPlotter import dataPlotter
from armDynamics import armDynamics
from ctrlPID import ctrlPID

# instantiate system, controller, and reference classes
arm = armDynamics()
controller = ctrlPID()
reference = signalGenerator(amplitude=20*np.pi/180.0, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = armAnimation()

t = P.t_start
y = arm.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = 0.5
        n = 0.0  #noise.random(t)
        u = controller.update(r, y + n)
        y = arm.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(arm.state)
    dataPlot.update(t, r, arm.state, u)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
