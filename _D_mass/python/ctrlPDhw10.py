import numpy as np
import massParam as P
from massDynamics import massDynamics

# set up the mass spring damper controller class
class ctrlPD:
    def __init__(self):
        #Givens
        zeta = .707 #damping ratio
        tr = 2.0 #rise time
        wn = 2.2/tr #omegan, natural frequency
        a0 = P.k/P.m
        a1 = P.b/P.m
        b0 = 1.0/P.m
        #PID gains
        self.kp = (wn**2 - a0) / b0
        self.kd = (2*zeta*wn - a1) / b0
        self.ki = 1.5 #integrator gain
        self.sigma = 0.05 #dirty derivative gain
        print('kp = ', self.kp)
        print('kd = ', self.kd)
        print('ki = ', self.ki)
        
        #integrator and differentiator
        self.zdot = 0.0 #estimated derivative of z
        self.z_d1 = 0.0 #z delayed by one sample
        self.error_dot = 0.0 #estimated derivative of error
        self.error_d1 = 0.0 #error delayed by one sample
        self.integrator = 0.0 #integrator
        
    def update(self, z_r, zCurrent): 
        # compute equilibrium force
        zerr = z_r - zCurrent
        #zdot = x[1][0]
        # integrate error
        self.integrator = self.integrator + (P.Ts/2.0)*(zerr + self.error_d1)
        # differentiate z using a discrete dirty derivative
        self.zdot = (2*self.sigma - P.Ts)/(2*self.sigma + P.Ts)*self.zdot \
            + (2.0/(2*self.sigma + P.Ts))*(zCurrent - self.z_d1)        
        F_tilde = self.kp*zerr + self.ki * self.integrator - self.kd*self.zdot
        Fout = self.saturate(F_tilde, P.F_max)
        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator = self.integrator + P.Ts/self.ki*(Fout - F_tilde)
        #Update the delayed variables
        self.error_d1 = zerr
        self.z_d1 = zCurrent
        return Fout
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u