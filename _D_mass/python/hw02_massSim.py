""" hw02_massSim.py
Jacob Child
Date Updated: Jan 19, 2023
HW Description: In the _D_mass directory, create the simulation
launch file hw02_massSim.py, and use it to animate the mass 
moving back and forth according to a sin wave with amplitude 
of 1 (m), frequency of 0.5 Hz, and y_offset of 0.2 (m).
"""
#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

#instantiate reference input classes
reference = signalGenerator(1.0, 0.5, 0.2) #amplitude, frequency, y_offset
thetaRef = signalGenerator(2.0*np.pi, 0.1) #amplitude, frequency
phiRef = signalGenerator(0.5, 0.1) #amplitude, frequency
tauRef = signalGenerator(5, 0.5) #amplitude, frequency

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end:
    # set variables
    r = reference.sin(t)
    theta = thetaRef.sin(t)
    phi = phiRef.sin(t)
    tau = tauRef.sawtooth(t)
    # update animation
    state = np.array([[theta], [phi], [0.0], [0.0]])
    animation.update(state)
    dataPlot.update(t, r, state, tau)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()