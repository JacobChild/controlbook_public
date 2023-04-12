import VTOLParam as P
from control import tf, bode
import control as cnt
import matplotlib.pyplot as plt
import numpy as np

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False

# Compute plant transfer functions
#compute the transfer function for F to H
FtoH = tf([1.0/(P.mc+2*P.mr)], [1.0, 0.0,0.0])
TautoTheta = tf([ 1.0/ (P.Jc+2*P.mr*P.d**2)], [1.0, 0.0, 0.0]) #? why was I short a 0? I think by accident
Fe = (P.mc + 2*P.mr)*P.g
#! wrong? ThetatoZ = tf([Fe/P.mu], [1.0, P.mu/(P.mc+2*P.mr), 0.0])
ThetatoZ = tf([Fe/(P.mc+2*P.mr)], [1.0, P.mu/(P.mc+2*P.mr), 0.0])
test1 = cnt.tf2ss(ThetatoZ)
#print(test1)

#plot the plant transfer functions
if __name__ == '__main__':
    # Bode plot for the plant
    fig = plt.figure()
    bode(FtoH, dB=dB_flag, margins=False)
    fig.axes[0].set_title('FtoH(s) for VTOL')
    fig = plt.figure()
    bode(TautoTheta, dB=dB_flag, margins=False)
    fig.axes[0].set_title('TautoTheta(s) for VTOL')
    fig = plt.figure()
    bode(ThetatoZ, dB=dB_flag, margins=False)
    fig.axes[0].set_title('ThetatoZ(s) for VTOL')

    # if you want specific values at specific frequencies, you can
    # do the following (but the magnitudes are absolute, not dB)
    mag, phase, omega = bode(FtoH, plot=False,
                             omega = [0.1, 10.0, 1000.0], dB=True)

    print('Close window to end program')
    plt.show()