""" hw02_massSim.py
Jacob Child
Date Updated: Jan 20, 2023
HW Description: In the _E_blockBeam directory, create the 
simulation launch file hw02_blockBeamSim.py, and use it to 
animate the block-on-beam system where the block is moving along
the beam according to a sin wave with amplitude of 0.05 (m), 
frequency of 0.5 Hz, and y_offset of 0.2 (m), and where the
angle of the beam is following a square wave with 
amplitude pi/8 (rad), and frequency of 0.1 Hz.
"""
#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter


#Instantiate reference input classes
reference = signalGenerator(.05, .5, .2) #amplitude, frequency, y_offset
blockPosRef = signalGenerator(.05, .5, .2) #amplitude, frequency, y_offset
beamAngleRef = signalGenerator(np.pi/8, .1) #amplitude, frequency
forceInputRef = None #this will be where the force input/controller is definded

#Instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end:
    # set variables
    r = reference.sin(t)
    blockPos = blockPosRef.sin(t)
    beamAngle = beamAngleRef.square(t)
    forceInput = forceInputRef #this will be where the force input/controller is definded
    # update animation
    state = np.array([[blockPos], [beamAngle], [0.0], [0.0]])
    animation.update(state)
    dataPlot.update(t, r, state, forceInput)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()