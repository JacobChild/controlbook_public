#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlStateFeedback import ctrlStateFeedback

#instantiate mass, controller, and reference classes
mass = massDynamics(alpha = 0.2) #Instantiated the massDynamics class as mass in this file
controller = ctrlStateFeedback() #Instantiated the ctrlStateFeedback class as controller in this file
# instantiate reference input classes
ZInputRef = signalGenerator(0.5, .08) #amplitude, frequency, y_offset

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
t = P.t_start #time starts at t_start
y = mass.h() #output of system at start of simulation
while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        zInput = ZInputRef.square(t)
        d = 0.0 #disturbance.step(t) #Get disturbance input
        n = 0.0 #noise input
        #update the dynamics
        zCurrent = y
        u  = controller.update(zInput, mass.state)
        y = mass.update(u + d) #propagate system
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, zInput, mass.state, u)
    #the pause causes the figure to be displayed during the simulation
    #plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()