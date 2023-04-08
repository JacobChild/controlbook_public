import massParam as P
from ctrlPIDhw10 import ctrlPID
import hw15bode as P15
from control import tf, bode_plot, bode
import matplotlib.pyplot as plt
import numpy as np

P10 = ctrlPID()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False # True

# Assign plant from previous homework solution
Plant = P15.Plant

# compute transfer function of controller
C_pid = tf([1.0/P.m *P10.kd, 1.0/P.m *P10.kp, 1.0/P.m * P10.ki], [1.0, P.b/P.m, P.k/P.m, 0.0])

#D16 part a
print("Tracking Error for unit ramp: ", (P.k)/(P10.ki))
    
if __name__ == '__main__':
    # display bode plots of transfer functions
    fig = plt.figure()
    bode_plot([Plant, Plant*C_pid], dB=dB_flag)
    fig.suptitle('Mass Spring Damper')
    plt.legend(['P(s)', 'C(s)P(s)'])
    #print('Close plot window to end program')
    #plt.show()
    
    
# D16 part b
#disturbence error percentage
magP, phaseP, omegaP = bode(Plant, plot=False,
                             omega = [0.1, 100.0], dB=dB_flag)
#for the controller
magC, phaseC, omegaC = bode(C_pid, plot=False,
                            omega = [0.1, 100.0], dB=dB_flag)
#convert to dB
magP = 20*np.log10(magP)
gammaDinPlant = 10**(-magP/20.0)
magC = 20*np.log10(magC)
gammaDinCtrl = 10**(-magC/20.0)

#difference
difference = magC - magP
#print("difference = ", difference)
gammaDin = 10**(-difference[0]/20.0)
print("gammaDin = ", gammaDin)
print("gammaN = ", 10**(magC[1]/20.0)) #? why not -magC[1]?