import matplotlib.pyplot as plt
import numpy as np
import armParam as P #has the arm Parameters
#Each of the above imports the given package/file as the given command
from signalGenerator import signalGenerator #imports the signalGenerator class
from armAnimation import armAnimation
from dataPlotter import dataPlotter
#each of the above imports a specific class/object from each of those files

# instantiate reference input classes
#instantiates (initiates) the class at each of the given values and saves them as different refs
reference = signalGenerator(amplitude=0.5, frequency=0.1)
thetaRef = signalGenerator(amplitude=2.0*np.pi, frequency=0.1)
tauRef = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = armAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    r = reference.square(t)
    theta = thetaRef.sin(t)
    tau = tauRef.sawtooth(t)
    # update animation
    state = np.array([[theta], [0.0]])  #state is made of theta, and theta_dot
    animation.update(state)
    dataPlot.update(t, r, state, tau)
    # advance time by t_plot
    t = t + P.t_plot  
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
