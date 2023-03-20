#Import all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from ctrlStateFeedbackwIntegrator import ctrlStateFeedback

#instantiate mass, controller, and reference classes
VTOL = VTOLDynamics(alpha = 0.2) #Instantiated the blockbeamDynamics class as blockbeam in this file
controller = ctrlStateFeedback()
#instantiate reference input classes
zPosRefSig = signalGenerator(3.0, 0.08, 2.5) #amplitude, frequency, y_offset
hInputRefSig = signalGenerator(5.0, 0.01, 5.0) #this will be where the torque input/controller is definded
disturbance = signalGenerator(amplitude=0.25)

#instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
t = P.t_start #time starts at t_start
y = VTOL.h()

while t < P.t_end: #main simulation loop
    #Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    #updates control and dynamics at faster simulation rate
    while t < t_next_plot:
        #Get referenced inputs from signal generators
        zPosRef = zPosRefSig.square(t)
        hPosRef = hInputRefSig.square(t)
        VTOLPosRef = np.array([[zPosRef], [hPosRef]])
        d = disturbance.step(t) #Get disturbance input
        n = 0.0
        x = y #this is the y from the past
        u = controller.update(VTOLPosRef, VTOL.state)
        #print(u)
        #update the dynamics
        y = VTOL.update(u+d)
        t = t + P.Ts #advance time by Ts
    #update animation and data plots
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, zPosRef,hPosRef,u[2][0], u[3][0]) #
    #the pause causes the figure to be displayed during the simulation
    #plt.pause(0.001)
    
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()