import matplotlib.pyplot as plt
import param as P
from dynamics import Dynamics
from controller import Controller
from signalGenerator import signalGenerator
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver

# instantiate arm, controller, and reference classes
plant = Dynamics(alpha=0.0)
controller = Controller()
reference = signalGenerator(amplitude=1.0, frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()

t = P.t_start  
y = plant.h()  
while t < P.t_end:  
    t_next_plot = t + P.t_plot
    while t < t_next_plot: 
        r = reference.square(t)
        d = 1.0   
        n = 0.  
        u, xhat, dhat = controller.update(r, y+n)  
        y = plant.update(u + d)  
        t = t + P.Ts  
    dataPlot.update(t, r, y, u)
    dataPlotObserver.update(t, plant.state, xhat, d, dhat)
    plt.pause(0.0001)  
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
