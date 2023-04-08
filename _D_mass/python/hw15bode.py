import numpy as np
import massParam as P
from control import tf, bode
import matplotlib.pyplot as plt

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False

# Compute plant transfer functions
Plant = tf([1.0/P.m],[1.0, P.b/P.m, P.k/P.m])
Ctrl = tf([1],[1.0, 0.0, 0.0])

if __name__ == '__main__':
    # Bode plot for the plant
    fig = plt.figure()
    bode(Plant, dB=dB_flag, margins=False)
    fig.axes[0].set_title('P(s) for mass')

    # if you want specific values at specific frequencies, you can
    # do the following (but the magnitudes are absolute, not dB)
    mag, phase, omega = bode(Plant, plot=False,
                             omega = [0.1, 10.0, 1000.0], dB=True)

    print('Close window to end program')
    plt.show()
    
magP, phase, omega = bode(Plant, plot=False,
                             omega = [0.1, 10.0, 1000.0])
#convert to dB
magP = 20*np.log10(magP)
print("mag = ", magP)