"""hw03_blockBeamSim.py
Jacob Child
Date Updated: Jan 27, 2023
HW Description: In the _E_blockBeam directory, create the file 
blockBeamDynamics.py and implement the dynamics equations of 
motion for the blockBeam system.  Use the file 
test_dynamics.py to test your dynamics file to ensure that it 
is providing the same outputs as my files for different 
initial conditions and inputs.  Create the file 
hw03_blockBeamSim.py, that simulates the blockBeam system 
when the input force is a sin wave with amplitude 0.5 (N), 
and frequency of 1 Hz, and a y_offset of 11.5 (N).  Note:  
The simulation will be unstable.
"""
#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockBeamDynamics import blockBeamDynamics

#instantiate mass, controller, and reference classes
blockbeam = blockBeamDynamics() #Instantiated the blockbeamDynamics class as blockbeam in this file
#instantiate reference input classes
forceInputRef = signalGenerator(0.5, 1.0, 11.5) #amplitude, frequency, y_offset
blockPosRefSig = signalGenerator(1.0, 0.1, 0.0) #amplitude, frequency, y_offset

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        forceInput = forceInputRef.sin(t)
        blockPosRef = blockPosRefSig.square(t)
        #update the dynamics
        y = blockbeam.update(forceInput)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockPosRef, blockbeam.state, forceInput) #! what do I put for the reference?
    #the pause causes the figure to be displayed during the simulation
    plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()