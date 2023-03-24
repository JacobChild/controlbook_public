import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import hummingbirdDynamics
from ctrlLatPID import ctrlLatPID

# instantiate pendulum, controller, and reference classes
hummingbird = hummingbirdDynamics(alpha=0.0)
controller = ctrlLatPID()
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.2)
psi_ref = SignalGenerator(amplitude=0.5, frequency=0.2)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[theta_ref.square(t)], [psi_ref.square(t)]])
        u, y_ref = controller.update(r, y)
        y = hummingbird.update(u)  # Propagate the dynamics
        t = t + P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    #force_plot = 
    #torque_plot =
    dataPlot.update(t, hummingbird.state, y_ref, u[0][0], u[1][0])

    # the pause causes figure to be displayed during simulation
    #plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
