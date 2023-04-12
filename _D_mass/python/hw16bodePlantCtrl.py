import massParam as P
from ctrlPIDhw10 import ctrlPID
import hw15bode as P15
from control import tf, bode_plot, bode
import matplotlib.pyplot as plt
import numpy as np

P10 = ctrlPID()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag # True

# Assign plant from previous homework solution
Plant = P15.Plant

# compute transfer function of controller
#! why is this wrong? C_pid = tf([1.0/P.m *P10.kd, 1.0/P.m *P10.kp, 1.0/P.m * P10.ki], [1.0, P.b/P.m, P.k/P.m, 0.0])
C_pid = tf([(P10.kd + P10.kp*P10.sigma), (P10.kp + P10.ki*P10.sigma), P10.ki], [P10.sigma, 1.0, 0.0])
#? it has something to do with the dirty derivative and sigmaS + 1 ?

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
magCP, phaseCP, omegaCP = bode(C_pid*Plant, plot=False,
                            omega = [0.1, 100.0], dB=dB_flag)
#if and elif statements to find gammaDin and gammaN if db_flag is true or false
if dB_flag == True:
    magP = 20*np.log10(magP)
    magCP = 20*np.log10(magCP)
    print("gammaDin = ", 10**((magP[0]-magCP[0]))/20)
    print("gammaN = ", 10**(magCP[1]/20))
elif dB_flag == False:
    print("gammaDin = ", magP[0]/magCP[0])
    print("gammaN = ", magCP[1])