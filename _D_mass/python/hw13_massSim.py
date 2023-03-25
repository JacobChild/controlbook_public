#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from massDynamics import massDynamics
from ctrlObserver import ctrlObserver
from dataPlotterObserver import dataPlotterObserver

#instantiate mass, controller, and reference classes
mass = massDynamics(alpha = 0.0) #Instantiated the massDynamics class as mass in this file
controller = ctrlObserver() #Instantiated 
# instantiate reference input classes
ZInputRef = signalGenerator(0.5, .04) #amplitude, frequency, y_offset
disturbance = signalGenerator(amplitude=0.25)

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
dataPlotObserver = dataPlotterObserver()
t = P.t_start #time starts at t_start
y = mass.h() #output of system at start of simulation
while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        zInput = ZInputRef.square(t)
        d = disturbance.step(t) #Get disturbance input
        n = 0.0 #noise.random(t) #noise input 
        #update the dynamics
        zCurrent = y
        u, xhat  = controller.update(zInput, zCurrent) 
        y = mass.update(u + d) #propagate system
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, zInput, mass.state, u)
    dataPlotObserver.update(t, mass.state, xhat)
    #the pause causes the figure to be displayed during the simulation
    #plt.pause(0.01)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()