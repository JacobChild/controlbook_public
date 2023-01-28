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
forceInputRef = signalGenerator(10.0, 1.0) #amplitude, frequency, y_offset

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        forceInput = forceInputRef.sin(t)
        #update the dynamics
        y = mass.update(forceInput)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, None, mass.state, forceInput)
    #the pause causes the figure to be displayed during the simulation
    plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()