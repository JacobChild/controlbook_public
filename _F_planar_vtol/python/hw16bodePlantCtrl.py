import numpy as np
from ctrlPIDhw10 import ctrlPD
import VTOLParam as P
from control import tf, bode, bode_plot
import matplotlib.pyplot as plt
import numpy as np
import hw15bode as P15

P10 = ctrlPD()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = P15.dB_flag

# Assign plants from previous homework solution
FtoH = P15.FtoH
TautoTheta = P15.TautoTheta
ThetatoZ = P15.ThetatoZ

# compute transfer function of the controllers
C_pid_FtoH = tf([(P10.kdh + P10.kph * P10.sigma), (P10.kph + P10.kih*P10.sigma), P10.kih], [P10.sigma,1.0,0.0])

C_pid_TautoTheta = tf([(P10.kdth + P10.kpth * P10.sigma), P10.kpth], [P10.sigma,1.0])

C_pid_ThetatoZ = tf([(P10.kdz + P10.kpz * P10.sigma), P10.kpz], [P10.sigma,1.0])

# plot them all with the plant on the same plot
if __name__ == '__main__':
    # display bode plots of transfer functions for the plants and controllers
    fig = plt.figure()
    bode_plot([FtoH, C_pid_FtoH], dB=dB_flag, margins=False)
    fig.suptitle('FtoH(s) for VTOL')
    plt.legend(['P(s)', 'C(s)P(s)'])
    fig = plt.figure()
    bode_plot([TautoTheta, C_pid_TautoTheta], dB=dB_flag, margins=False)
    fig.suptitle('TautoTheta(s) for VTOL')
    plt.legend(['P(s)', 'C(s)P(s)'])
    fig = plt.figure()
    bode_plot([ThetatoZ, C_pid_ThetatoZ], dB=dB_flag, margins=False)
    fig.suptitle('ThetatoZ(s) for VTOL')
    plt.legend(['P(s)', 'C(s)P(s)'])
    #plt.show()
    

#Problem F16 part b FtoH
mag, phase, omega = bode(FtoH*C_pid_FtoH, plot=False,
                             omega = [30.0], dB = dB_flag)
if dB_flag == True:
    mag = 20*np.log10(mag)
    print("gammaN_FtoH = ", 10.0**(mag/20.0))
elif dB_flag == False:
    print("gammaN_FtoH = ", mag)

#part c TautoTheta input disturbance, so difference between plant
magCP, phase, omega = bode(C_pid_TautoTheta*TautoTheta, plot=False,
                         omega = [2.0], dB = dB_flag)
magP, phaseP, omegaP = bode(TautoTheta, plot=False,
                            omega = [2.0], dB = dB_flag)
if dB_flag == True:
    mag = 20*np.log10(magCP)
    magP = 20*np.log10(magP)
    gammaDin_TautoTheta = 10.0**((magP - magCP)/20.0)
    print("gammaDin_TautoTheta = ", gammaDin_TautoTheta)
elif dB_flag == False:
    gammaDin_TautoTheta = magP/magCP
    print("gammaDin_TautoTheta = ", gammaDin_TautoTheta)

#part d find the sensor range for theta with measurement noise less than .1 deg
rad = .1#*np.pi/180.0
term1 = 1.0/(P.Jc + 2.0*P.mr*P.d**2)
#print("term1 = ", term1)
root1 = (-P10.kdth + np.sqrt(P10.kdth**2 - 4*-P10.kpth*rad*(P.Jc + 2.0*P.mr*P.d**2)))/(2*rad*(P.Jc + 2.0*P.mr*P.d**2))
root2 = (-P10.kdth - np.sqrt(P10.kdth**2 - 4*-P10.kpth*rad*(P.Jc + 2.0*P.mr*P.d**2)))/(2*rad*(P.Jc + 2.0*P.mr*P.d**2))
print("bounds for theta = ", root1, "and ", root2)

#part e ThetatoZ tracking error gammar
magCP, phase, omega = bode(C_pid_ThetatoZ*ThetatoZ, plot=False,
                         omega = [0.1, 0.01], dB = dB_flag)
magP, phaseP, omegaP = bode(ThetatoZ, plot=False,
                            omega = [0.1, 0.01], dB = dB_flag)
if dB_flag == True:
    magCP = 20*np.log10(magCP)
    print("gammaDout_ThetatoZ = ", 1.0/magCP[0])
    print("gammaRef_ThetatoZ = ", 1.0/magCP[1])
elif dB_flag == False:
    print("gammaDout_ThetatoZ = ", 1.0/magCP[0])
    print("gammaRef_ThetatoZ = ", 1.0/magCP[1])

