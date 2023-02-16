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
    u = np.array([[5.0],
                  [5.0]])
    #update the dynamics
    y = hummingbird.update(u)
    state = hummingbird.state
    ref = np.array([[phi], [theta], [psi]])
    animation.update(t, state)
    pwm_left = u[0][0]
    pwm_right = u[1][0]
    force = P.km * (pwm_left + pwm_right)
    torque = P.d * P.km * (pwm_left - pwm_right)
    dataPlot.update(t, state, ref, force, torque)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()