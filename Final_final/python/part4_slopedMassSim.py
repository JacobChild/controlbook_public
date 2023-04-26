import matplotlib.pyplot as plt
import numpy as np
import slopedMassParam as P
from signalGenerator import signalGenerator
from slopedMassAnimation import slopedMassAnimation
from dataPlotter import dataPlotter
from slopedMassDynamics import slopedMassDynamics
from ctrlLoop import ctrlLoop

# instantiate system, controller, and reference classes
slopedMass = slopedMassDynamics()
controller = ctrlLoop('digital_filter')
reference = signalGenerator(amplitude=0.5, frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = slopedMassAnimation()

t = P.t_start
y = slopedMass.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = 1.0
        n = 0.0  #noise.random(t)
        u = controller.update(r, y + n)
        y = slopedMass.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(slopedMass.state)
    dataPlot.update(t, r, slopedMass.state, u)
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
