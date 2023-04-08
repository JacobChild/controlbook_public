import numpy as np
from ctrlPIDhw10 import ctrlPD
import VTOLParam as P
from control import tf, bode, bode_plot
import matplotlib.pyplot as plt
import numpy as np
import hw15bode as P15

P10 = ctrlPD()

# flag to define if using dB or absolute scale for M(omega)
dB_flag = False

# Assign plants from previous homework solution
FtoH = P15.FtoH
TautoTheta = P15.TautoTheta
ThetatoZ = P15.ThetatoZ

# compute transfer function of the controllers
term = 1.0/(P.mc + 2*P.mr)
C_pid_FtoH = tf([term*P10.kdh, term*P10.kph, term*P10.kih], [1.0, 0.0, 0.0])
term2 = 1.0/(P.Jc + 2*P.mr*P.d**2)
C_pid_TautoTheta = tf([term2*P10.kdth, term2*P10.kpth], [1.0, 0.0])
term3 = 1.0/(P.mc + 2*P.mr)
C_pid_ThetatoZ = tf([-P15.Fe*term3*P10.kdz, -P15.Fe*term3*P10.kpz, -P15.Fe*term3*P10.kiz], [1.0, P.mu*term3, 0.0])

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
mag, phase, omega = bode(C_pid_FtoH, plot=False,
                             omega = [30.0], dB = dB_flag)
gammaN_FtoH = mag
print('gammaN_FtoH = ', gammaN_FtoH)

#part c TautoTheta input disturbance, so difference between plant
mag, phase, omega = bode(C_pid_TautoTheta, plot=False,
                         omega = [2.0], dB = dB_flag)
magP, phaseP, omegaP = bode(TautoTheta, plot=False,
                            omega = [2.0], dB = dB_flag)
#convert to dB
mag = 20*np.log10(mag)
magP = 20*np.log10(magP)
gammaDin_TautoTheta = 10**(-(mag - magP)/20.0)
print("gammaDin_TautoTheta = ", gammaDin_TautoTheta)

#part d find the sensor range for theta with measurement noise less than .1 deg
rad = .1#*np.pi/180.0
term1 = 1.0/(P.Jc + 2.0*P.mr*P.d**2)
#print("term1 = ", term1)
root1 = (-P10.kdth + np.sqrt(P10.kdth**2 - 4*-P10.kpth*rad*(P.Jc + 2.0*P.mr*P.d**2)))/(2*rad*(P.Jc + 2.0*P.mr*P.d**2))
root2 = (-P10.kdth - np.sqrt(P10.kdth**2 - 4*-P10.kpth*rad*(P.Jc + 2.0*P.mr*P.d**2)))/(2*rad*(P.Jc + 2.0*P.mr*P.d**2))
print("bounds for theta = ", root1, "and ", root2)

#part e ThetatoZ tracking error gammar
mag, phase, omega = bode(C_pid_ThetatoZ, plot=False,
                         omega = [0.1, 0.01], dB = dB_flag)
magP, phaseP, omegaP = bode(ThetatoZ, plot=False,
                            omega = [0.1, 0.01], dB = dB_flag)
gammaRef_ThetatoZ = 1.0 / mag[0]
print("gammaRef_ThetatoZ = ", gammaRef_ThetatoZ)
gammaDout_ThetatoZ = 1.0 / mag[1]
print("gammaDout_ThetatoZ = ", gammaDout_ThetatoZ)
