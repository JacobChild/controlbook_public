#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics
from ctrlStateFeedback import ctrlStateFeedback

#instantiate mass, controller, and reference classes
blockbeam = blockBeamDynamics(alpha = 0.2) #Instantiated the blockbeamDynamics class as blockbeam in this file
controller = ctrlStateFeedback()
#instantiate reference input classes
#ForceInputRef = signalGenerator(0.5, 1.0, 11.5) #amplitude, frequency, y_offset
#blockPosRefSig = signalGenerator(.125, 0.02, 0.250) #amplitude, frequency, y_offset
blockPosRefSig = signalGenerator(.125, 0.08, 0.250)

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
        n = 0.0
        xCurrent = y #in the form [z][theta], this y is from the past
        u = controller.update(blockPosRef, blockbeam.state)
        #update the dynamics
        y = blockbeam.update(u)
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