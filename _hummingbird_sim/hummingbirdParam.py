# Hummingbird Parameter File
import numpy as np
# Initial Conditions
phi0 = 0.0 * np.pi / 180  # roll angle in rads
theta0 = 0 * np.pi / 180  # pitch angle in rads
psi0 = 0.0 * np.pi / 180  # yaw angle in rads
phidot0 = 0.0              # roll rate in rads/sec
thetadot0 = 0.0         # pitch rate in rads/sec
psidot0 = 0.0              # yaw rate in rads/sec
# Physical parameters of the hummingbird known to the controller
g = 9.81 # m/s^2
ell1 = .2475 # m
ell2 = -.039 # m
ell3x = -.007 # m
ell3y = -.007 # m
ell3z = .018 # m
ellT = .355
d = .12
m1 = .108862
J1x = .000189 # kg*m^2
J1y = .001953
J1z = .001894
m2 = .4717
J2x = .00231 # kg*m^2
J2y = .003274
J2z = .003416
m3 = .1905
J3x = .0002222
J3y = .0001956
J3z = .000027
km = g * (m1 * ell1 + m2 * ell2) / ellT  # need to find this experimentally for hardware


# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 30.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.1  # the plotting and animation is updated at this rate
# saturation limits
force_max = 2.0                # Max force N
torque_max = 5.0                # Max torque, Nm

