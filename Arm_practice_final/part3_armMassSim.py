import matplotlib.pyplot as plt
import numpy as np
import armParam as P
from signalGenerator import signalGenerator
from armAnimation import armAnimation
from dataPlotter import dataPlotter
from armDynamics import armDynamics
from ctrlObsv import ctrlObsv
from dataPlotterObserver import dataPlotterObserver

# instantiate system, controller, and reference classes
arm = armDynamics()
controller = ctrlObsv()
reference = signalGenerator(amplitude=20*np.pi/180.0, frequency=0.1)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = armAnimation()

t = P.t_start
y = arm.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = 0.5
        n = 0.0  #noise.random(t)
        u, xhat, dhat = controller.update(r, y + n)
        y = arm.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(arm.state)
    dataPlot.update(t, r, arm.state, u)
    dataPlotObserver.update(t, arm.state, xhat, d, dhat)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
