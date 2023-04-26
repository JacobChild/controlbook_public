import matplotlib.pyplot as plt
import numpy as np
import slopedMassParam as P
from signalGenerator import signalGenerator
from slopedMassAnimation import slopedMassAnimation
from dataPlotter import dataPlotter
from slopedMassDynamics import slopedMassDynamics
from ctrlObsv import ctrlObsv
from dataPlotterObserver import dataPlotterObserver

# instantiate system, controller, and reference classes
slopedMass = slopedMassDynamics()
controller = ctrlObsv()
reference = signalGenerator(amplitude=0.5, frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = slopedMassAnimation()

t = P.t_start
y = slopedMass.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = 1.0
        n = 0.0  #noise.random(t)
        u, xhat, dhat = controller.update(r, y + n)
        y = slopedMass.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(slopedMass.state)
    dataPlot.update(t, r, slopedMass.state, u)
    dataPlotObserver.update(t, slopedMass.state, xhat, d, dhat)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
