""" hw03_VTOLSim.py
Jacob Child
Date Updated: Jan 27, 2023
HW Description: In the _F_VTOL directory, create the file 
VTOLDynamics.py and implement the dynamics equations of 
motion for the planar VTOL system.  Use the file 
test_dynamics.py to test your dynamics file to ensure that 
it is providing the same outputs as my files for different 
initial conditions and inputs.  Create the file 
hw03_VTOLSim.py, that simulates the planar VTOL system when 
the thrust on both motors is a sin wave with amplitude 10 (N), 
and frequency of 1 Hz.  Note:  The simulation will not do 
anything interesting.
"""

#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics

#instantiate mass, controller, and reference classes
VTOL = VTOLDynamics() #Instantiated the VTOLDynamics class as VTOL in this file
#instantiate reference input classes
forceRightInputRef = signalGenerator(10.0, 1.0, 0.0) #amplitude, frequency, y_offset
forceLeftInputRef = signalGenerator(10.0, 1.0, 0.0) #amplitude, frequency, y_offset
VTOLPosRefSig = signalGenerator(1.0, 0.1, 0.0) #amplitude, frequency, y_offset

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        forceRightInput = forceRightInputRef.sin(t)
        forceLeftInput = forceLeftInputRef.sin(t)
        VTOLPosRef = VTOLPosRefSig.square(t)
        #update the dynamics
        forceInput = np.array([[forceRightInput], [forceLeftInput]])
        y = VTOL.update(forceInput)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(VTOL.state)
    dataPlot.update(t, VTOLPosRef, VTOL.state, forceInput) 
    # the pause causes the figure to be displayed during the simulation
    plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
    