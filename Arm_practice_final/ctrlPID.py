import numpy as np
import armParam as P


class ctrlPID:
    def __init__(self):
        tr = 0.25 #sec
        wn = 2.2/tr
        zeta = 0.707
        a1 = 0.0
        b0 = 3.0 / (P.m * P.ell**2)
        a0 = 0.0
        self.kd = (2.0*zeta*wn - a1) / b0 #these are general equations and should work for all PD systems
        self.kp = (wn**2 - a0) / b0 
        self.ki = 2.0 #Integrator gain that I tune
        print("kd: ", self.kd, " kp: ", self.kp, " ki: ", self.ki)
        #other needed parameters
        self.sigma = 0.005
        self.Ts = P.Ts
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.limit = P.tau_max #his built in saturation function uses self.limit
        #variables and delayed variables for calculation
        self.thetadot = 0.0
        self.integrator = 0.0
        self.error_d1 = 0.0
        self.theta_d1 = 0.0 #delayed theta

    def update(self, theta_r, y):
        theta = y
        error = theta_r - theta
        #integrate on error
        self.integrator = self.integrator + (P.Ts/2.0)*(error + self.error_d1)
        #compute derivative
        self.thetadot = self.beta*self.thetadot + (1.0-self.beta) * ((theta - self.theta_d1) / P.Ts)
        tau_tilde = self.kp * error - self.kd * self.thetadot + self.ki * self.integrator
        th_eq = 0.0
        tau_eq = P.m * P.g * P.ell / 2.0 * np.cos(th_eq)
        tau = self.saturate(tau_tilde + tau_eq)
        #integrator anti windup
        if self.ki != 0.0:
            self.integrator =  self.integrator + P.Ts/self.ki*(tau - (tau_tilde+tau_eq)) #?ie if it is saturating decrease the integrator
        
        #update delayed variables
        self.error_d1 = error
        self.theta_d1 = theta
        
        
        return tau.item(0)

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u







