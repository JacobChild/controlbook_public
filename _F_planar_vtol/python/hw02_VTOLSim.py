""" hw02_VTOLSim.py
Jacob Child
Date Updated: Jan 20, 2023
HW Description: In the _F_VTOL directory, create the simulation launch file 
hw02_VTOLSim.py, and use it to animate the planar VTOL system where the the 
position z follows a sin wave with amplitude of 4 (m), frequency of 0.1 Hz, 
and y_offset of 5 (m), altitude h follows a square wave with amplitude 
of 2 (m), frequency of 0.1 Hz,  and y_offset of 2 (m), and the 
angle theta follows a sin wave with amplitude pi/8 (rad), and 
frequency of 0.5 Hz
"""

#! Look at the animation & plotter files
#IMport all of the needed classes and files
import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter

#Instantiate reference input classes
VTOLZPosAct = signalGenerator(4, .1, 5) #amplitude, frequency, y_offset
VTOLZPosRef = signalGenerator(4, .1, 5) #amplitude, frequency, y_offset
AltitudeAct = signalGenerator(2, .1, 2.0) #amplitude, frequency, y_offset
AltitudeRef = signalGenerator(2, .1, 2.0) #amplitude, frequency 
thetaRef = signalGenerator(np.pi/8, .5) #amplitude, frequency
forceInputRef = 0.0 #this will be where the force input/controller is definded
torqueInputRef = 0.0 #this will be where the torque input/controller is definded

#Instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()
t = P.t_start #time starts at t_start
while t < P.t_end:
    # set variables
    VTOLZPos = VTOLZPosAct.sin(t)
    Altitude = AltitudeAct.square(t)
    theta = thetaRef.sin(t)
    VTOLZCurRef = VTOLZPosRef.square(t)
    AltitudeCurRef = AltitudeRef.square(t)
    forceInput = forceInputRef #this will be where the force input/controller is definded
    torqueInput = torqueInputRef #this will be where the torque input/controller is definded
    # update animation
    state = np.array([[VTOLZPos], [Altitude], [theta], [0.0]])
    animation.update(state, AltitudeCurRef)
    dataPlot.update(t, state, VTOLZCurRef, AltitudeCurRef, forceInput, torqueInput)
    # advance time by t_plot
    t = t + P.t_plot
    plt.pause(0.001)
   
# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()