#us this file to show equilibrium ie it is ctrlequilibrium.py

import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import hummingbirdDynamics

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)
force_ref = SignalGenerator(amplitude=0.5, frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingbird = hummingbirdDynamics()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    phi = phi_ref.sin(t)
    theta = theta_ref.sin(t)
    psi = psi_ref.sin(t)
    Fe = (P.m1*P.ell1 + P.m2*P.ell2)*P.g/ (P.ellT)
    F = Fe
    tau = 0.0
    Km = P.km
    ul = 1.0/(2.0*Km)*(F + tau/P.d)
    ur = 1.0/(2.0*Km)*(F - tau/P.d)
    u = np.array([[ul],
                  [ur]])
    #update the dynamics
    y = hummingbird.update(u)
    state = hummingbird.state
    ref = np.array([[phi], [theta], [psi]])
    animation.update(t, state)
    force = P.km * ul + P.km * ur
    torque = P.d * P.km * (P.d *ul - P.d * ur)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()