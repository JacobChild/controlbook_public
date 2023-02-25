import matplotlib.pyplot as plt
import numpy as np
import param as P
from signalGenerator import signalGenerator
from dataPlotter import dataPlotter
from dynamics import Dynamics
from ctrlPD import ctrlPD

# instantiate system, controller, and reference classes
system = Dynamics()
controller = ctrlPD()
reference = signalGenerator(amplitude=5, frequency=0.1)
dataPlot = dataPlotter()

t = P.t_start  
y = system.h()  
while t < P.t_end:  
    t_next_plot = t + P.t_plot
    while t < t_next_plot: 
        r = reference.square(t)
        x = system.state
        u = controller.update(r, x)  
        y = system.update(u)  
        t = t + P.Ts 
    dataPlot.update(t, r, system.state, u)
    plt.pause(0.0001)  
print('kp: ', controller.kp)
print('kd: ', controller.kd)        
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
