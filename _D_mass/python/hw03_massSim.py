"""hw03_massSim.py
Jacob Child
Date Updated: Jan 26, 2023
HW Description: In the _D_mass directory, create the file 
massDynamics.py and implement the dynamics equations of motion
for the mass-spring-damper system.  Use the file 
test_dynamics.py to test your dynamics file to ensure that 
it is providing the same outputs as my files for different 
initial conditions and inputs.  Create the file 
hw03_massSim.py, that simulates the mass system when the 
input force is a sin wave with amplitude 10 (N), and 
frequency of 1 Hz.  Note:  The simulation will not do 
anything interesting.
"""
#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics

#instantiate mass, controller, and reference classes
mass = massDynamics() #Instantiated the massDynamics class as mass in this file
#instantiate reference input classes
reference = signalGenerator(0.0, 0.5, 0.2) #amplitude, frequency, y_offset
massPosRef = signalGenerator(1.0, .5, .2) #amplitude, frequency, y_offset
forceInputRef = signalGenerator(10.0, 1.0, 0.0) #amplitude, frequency, y_offset

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end: #! compare this while loop to the one in hw03_armSim.py
    #set variables
    r = reference.sin(t)
    massPos = massPosRef.sin(t)
    forceInput = forceInputRef.sin(t)
    #update animation
    state = np.array([[massPos], [None], [None], [None]])
    animation.update(state)
    dataPlot.update(t, r, state, forceInput)
    #advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()