import numpy as np
import slopedMassParam as P


class ctrlPID:
    def __init__(self):
        tr = 0.5 #sec
        wn = 2.2/tr
        zeta = 2.0
        a1 = P.b/ (P.m)
        b0 = 1.0 / (P.m)
        a0 = P.k1 / (P.m)
        self.kd = (2.0*zeta*wn - a1) / b0 #these are general equations and should work for all PD systems
        self.kp = (wn**2 - a0) / b0 
        self.ki = 5.0 #Integrator gain that I tune
        print("kd: ", self.kd, " kp: ", self.kp, " ki: ", self.ki)
        #other needed parameters
        self.sigma = P.sigma #0.05 I believe
        self.Ts = P.Ts
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.limit = P.F_max #his built in saturation function uses self.limit
        #variables and delayed variables for calculation
        self.zdot = 0.0
        self.integrator = 0.0
        self.error_d1 = 0.0
        self.z_d1 = 0.0 #delayed z

    def update(self, z_r, z):
        z = z #[0][0]
        error = z_r - z
        #integrate on error
        self.integrator = self.integrator + (P.Ts/2.0)*(error + self.error_d1)
        #compute derivative
        self.zdot = self.beta*self.zdot + (1.0-self.beta) * ((z - self.z_d1) / P.Ts)
        F_tilde = self.kp * error - self.kd * self.zdot + self.ki * self.integrator
        
        z_eq = 0.0 
        F_eq = P.k1 * z_eq + P.k2 * z_eq**3 - P.m * P.g * np.sin(np.pi/4)
        F = self.saturate(F_tilde+F_eq)
        #integrator anti windup
        if self.ki != 0.0:
            self.integrator =  self.integrator + P.Ts/self.ki*(F - (F_tilde+F_eq)) #?ie if it is saturating decrease the integrator
        
        #update delayed variables
        self.error_d1 = error
        self.z_d1 = z
        return F
    
    
    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u



    #this is the saturate function he gave us, I will use the one from the practice final
    # def saturate(u, limit):
    #     if abs(u) > limit:
    #         u = limit*np.sign(u)
    #     return u







