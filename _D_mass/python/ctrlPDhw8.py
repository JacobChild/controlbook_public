import numpy as np
import massParam as P
from massDynamics import massDynamics

# set up the mass spring damper controller class
class ctrlPD:
    def __init__(self):
        #Givens
        zeta = .707 #damping ratio
        tr = 0.50 #rise time
        wn = 2.2/tr #omegan, natural frequency
        a0 = P.k/P.m
        a1 = P.b/P.m
        b0 = 1.0/P.m
        #PD gains
        self.kp = (wn**2 - a0) / b0
        self.kd = (2*zeta*wn - a1) / b0
        print('kp = ', self.kp)
        print('kd = ', self.kd)
        
    def update(self, z_r, x): 
        # compute equilibrium force
        zerr = z_r - x[0][0]
        zdot = x[1][0]
        #Calculate F equilibrium
        #Fe = P.k * x[0][0]
        # compute the linearized force using PD control
        F_tilde = self.kp*zerr - self.kd*zdot
       # FNew = massDynamics.saturate(F_tilde, P.F_max)
        # compute total force
        Fout = F_tilde #+ Fe
        return Fout