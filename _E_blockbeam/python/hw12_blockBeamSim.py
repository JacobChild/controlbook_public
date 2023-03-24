#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics
from ctrlStateFeedbackwIntegrator import ctrlStateFeedbackwIntegrator

#instantiate mass, controller, and reference classes
blockbeam = blockBeamDynamics(alpha = 0.0) #Instantiated the blockbeamDynamics class as blockbeam in this file
controller = ctrlStateFeedbackwIntegrator()
blockPosRefSig = signalGenerator(amplitude=.125, frequency=0.05, y_offset = 0.25)
disturbance = signalGenerator(amplitude=0.25, frequency = 0.0)

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

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
        n = 0.0
        #xCurrent = blockbeam.state #in the form [z][theta], this y is from the past
        u = controller.update(blockPosRef, blockbeam.state)
        #update the dynamics
        y = blockbeam.update(u+d)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockPosRef, blockbeam.state, u) #
    #the pause causes the figure to be displayed during the simulation
    #plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()