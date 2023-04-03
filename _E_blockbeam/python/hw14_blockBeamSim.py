#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics
from ctrlDisturbanceObserver import ctrlDisturbanceObserver
from dataPlotterObserver import dataPlotterObserver

#instantiate mass, controller, and reference classes
blockbeam = blockBeamDynamics(alpha = 0.2) #Instantiated the blockbeamDynamics class as blockbeam in this file
controller = ctrlDisturbanceObserver()
blockPosRefSig = signalGenerator(amplitude=.125, frequency=0.05, y_offset = 0.25)
disturbance = signalGenerator(amplitude=0.25, frequency = 0.0)
noise_z = signalGenerator(amplitude=0.00)
noise_th = signalGenerator(amplitude=0.00)

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
dataPlotObserver = dataPlotterObserver()

t = P.t_start #time starts at t_start
y = blockbeam.h()

while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        #forceInput = forceInputRef.sin(t)
        blockPosRef = blockPosRefSig.square(t)
        d = disturbance.step(t) #Get disturbance input
        #n = np.array([[0.0], [0.0]])
        n = np.array([[noise_z.random(t)], [noise_th.random(t)]])
        #xCurrent = blockbeam.state #in the form [z][theta], this y is from the past
        u, xhat, dhat = controller.update(blockPosRef, y + n)
        #update the dynamics
        y = blockbeam.update(u+d)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockPosRef, blockbeam.state, u) 
    dataPlotObserver.update(t, blockbeam.state, xhat, d, dhat)
    #the pause causes the figure to be displayed during the simulation
    #plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()