# Inverted Pendulum Parameter File
import numpy as np
# import control as cnt

# Physical parameters of the arm known to the controller
m = 0.5     # Mass of the arm, kg
ell = 0.3    # Length of the arm, m
g = 9.8       # Gravity, m/s**2
b = 0.01      # Damping coefficient, Nms

# parameters for animation
length = 1.    # length of arm in animation
width = 0.3   # width of arm in animation

# Initial Conditions
theta0 = 0.0*np.pi/180  # ,rads
thetadot0 = 0.0         # ,rads/s

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 20.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate

# saturation limits
tau_max = 1.0                # Max torque, N-m

