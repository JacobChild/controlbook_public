import rodMassParam as P
import matplotlib.pyplot as plt
from control import TransferFunction as tf
from control import tf, bode, margin, step_response, mag2db
import numpy as np
import loopshape_tools as ls
from ctrlPID import ctrlPID
PID = ctrlPID()

# Compute plant transfer functions
Plant = tf([1.0/(P.m * P.ell**2)], #numerator
           [1.0, P.b/(P.m * P.ell**2), P.k1 / (P.m * P.ell**2)]) #this comes from the plant, make sure each term has something, even if a 0.0
C_pid = tf([(PID.kd+PID.kp*PID.sigma), 
            (PID.kp+PID.ki*PID.sigma), 
            PID.ki],
           [PID.sigma, 1, 0]) # this should be the same for every PID controller I believe
PLOT = True
#PLOT = False
dB_flag = True

#######################################################################
#   Control Design
#######################################################################
C = C_pid * ls.lead(w=22.149,M=10.0) * ls.lag(z=30.0, M =25.0) * ls.lpf(p=100.0)
#lead is for phase margin, lag is for disturbances/tracking, lpf is for noise

###########################################################
# add a prefilter to eliminate the overshoot
###########################################################
F = tf(1, 1) * ls.lpf(p=5.0) #originally this was p=3.0 from the example I think?
##############################################
#  Convert Controller to State Space Equations if following method in 18.1.7
##############################################
C_num = np.asarray(C.num[0])
C_den = np.asarray(C.den[0])
F_num = np.asarray(F.num[0])
F_den = np.asarray(F.den[0])

if __name__ == "__main__":
    # calculate bode plot and gain and phase margin for just the plant dynamics
        #for the above see the quick code just below
        
    #### Code added to find gammaN and gammaR and to plot the noise and tracking specifications
    #for the controller
    #Also quick code to plot just the plant
    mag, phase, omega = bode(Plant, dB=True,
                             omega=np.logspace(-3, 5),
                             Plot=True, label="$P(s)$")

    gm, pm, Wcg, Wcp = margin(Plant * C_pid)
    magCP, phaseCP, omegaCP = bode(Plant*C_pid, plot=False,
                            omega = [0.001, 100.0], dB=dB_flag) #TODO fill out these omega's for gammaN and gammaR
    mag4Plt, phase4Plt, omega4Plt = bode(Plant*C_pid, plot=False,
                            omega = [0.02, 2000.0], dB=dB_flag) #TODO fill out these omegas for the tracking and noise specifications
    
    #Tracking and noise specifications
    ls.spec_tracking(gamma=0.1*1.0/mag4Plt[0], omega=0.02, flag=dB_flag)
    ls.spec_noise(gamma=0.1*mag4Plt[1], omega=2000.0, flag=dB_flag)
   
    
    print("MagCP: ", magCP)
    print("for original C_pid system:")
    #It will spit out absolute magnitude, so I will not need to convert
    print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", gm, " Wcg: ", Wcg)
    print("gammaR = ", 1.0/magCP[0])
    print("gammaN = ", magCP[1])
    #A few additional print statements
    print("\n")
    print("C(s): ", C)
    print("\n")
    print("F(s): ", F)
    
    
    # calculate bode plot and gain and phase margin for original PID * plant dynamics
    mag, phase, omega = bode(Plant * C_pid, dB=True,
                             omega=np.logspace(-3, 5),
                             Plot=True, label="$P(s)C_{pid}(s)$")

    gm, pm, Wcg, Wcp = margin(Plant * C_pid)
    print("for original C_pid system:")
    if dB_flag is True:
        print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", mag2db(gm), " Wcg: ", Wcg)
    elif dB_flag is False:
        print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", gm, " Wcg: ", Wcg)

    #########################################
    #   Define Design Specifications
    #########################################

    # plot the effect of adding the new compensator terms
    mag, phase, omega = bode(Plant * C, dB=dB_flag,
                             omega=np.logspace(-4, 5),
                             plot=True, label="$P(s)C_{final}(s)$")

    gm, pm, Wcg, Wcp = margin(Plant * C)
    print("for final P*C:")
    if dB_flag is True:
        print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", mag2db(gm), " Wcg: ", Wcg)
    elif dB_flag is False:
        print(" pm: ", pm, " Wcp: ", Wcp, "gm: ", gm, " Wcg: ", Wcg)

    plt.figure(1)
    fig = plt.gcf()
    fig.axes[0].legend()

    ############################################
    # now check the closed-loop response with prefilter
    ############################################
    # Closed loop transfer function from R to Y - no prefilter
    CLOSED_R_to_Y = (Plant * C / (1.0 + Plant * C))
    # Closed loop transfer function from R to Y - with prefilter
    CLOSED_R_to_Y_with_F = (F * Plant * C / (1.0 + Plant * C))
    # Closed loop transfer function from R to U - no prefilter
    CLOSED_R_to_U = (C / (1.0 + Plant * C))
    # Closed loop transfer function from R to U - with prefilter
    CLOSED_R_to_U_with_F = (F * C / (1.0 + Plant * C))

    fig = plt.figure(2)
    plt.clf()
    plt.grid(True)
    mag, phase, omega = bode(CLOSED_R_to_Y, dB=dB_flag, plot=True,
                             color=[0, 0, 1], label='closed-loop $\\frac{Y}{R}$ - no pre-filter')
    mag, phase, omega = bode(CLOSED_R_to_Y_with_F, dB=dB_flag, plot=True,
                             color=[0, 1, 0], label='closed-loop $\\frac{Y}{R}$ - with pre-filter')
    fig.axes[0].set_title('Closed-Loop Bode Plot')
    fig.axes[0].legend()

    plt.figure(4)
    plt.clf()
    plt.subplot(211), plt.grid(True)
    T = np.linspace(0, 2, 100)
    _, yout_no_F = step_response(CLOSED_R_to_Y, T)
    _, yout_F = step_response(CLOSED_R_to_Y_with_F, T)
    plt.plot(T, yout_no_F, 'b', label='response without prefilter')
    plt.plot(T, yout_F, 'g', label='response with prefilter')
    plt.legend()
    plt.ylabel('Step Response')

    plt.subplot(212), plt.grid(True)
    _, Uout_no_F = step_response(CLOSED_R_to_U, T)
    _, Uout_F = step_response(CLOSED_R_to_U_with_F, T)
    plt.plot(T, Uout_no_F, color='b', label='control effort without prefilter')
    plt.plot(T, Uout_F, color='g', label='control effort with prefilter')
    plt.ylabel('Control Effort')
    plt.legend()

    plt.show()

