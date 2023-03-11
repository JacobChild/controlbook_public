import numpy as np
import blockbeamParam as P

# set up the block beam controller class
class ctrlPD:
    def __init__(self):
        #Givens
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_z = 2.0 # rise time for outer loop
        zeta_th = 0.707  # inner loop damping ratio 
        M = 10.0  # Time scale separation between loops
        tr_th = tr_z / M  # rise time for inner loop
        zeta_z = 0.707  # outer loop damping ratio
        self.sigma = 0.05  # cutoff freq for dirty derivative
        # saturation limits
        self.force_max = P.Fmax  
            # maximum commanded base angle
        #---------------------------------------------------
        #                    Inner Loop
        #---------------------------------------------------
        # PD givens for inner loop
        wn_th = 2.2 / tr_th
        a0_th = 0.0
        a1_th = 0.0
        z0_err_th = P.length/2.0
        b0_th = P.length/(P.m2*P.length**2 / 3.0 + P.m1 * z0_err_th**2)
        #PID gains
        self.kpth = (wn_th**2 - a0_th) / b0_th
        self.kdth = (2.0*zeta_th*wn_th - a1_th) / b0_th
        print('kpth = ', self.kpth)
        print('kdth = ', self.kdth)
        
        # inner loop integrator and differentiator values
        self.thetadot = 0.0 #estimated derivative of theta
        self.theta_d1 = 0.0 #theta delayed by one sample
        self.error_thetadot = 0.0 #estimated derivative of error
        self.integrator_theta = 0.0 #integrator
        
        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # PD design for outer loop
        wn_z =2.2 / tr_z
        a0_z = 0.0
        a1_z = 0.0
        z0_err_z = 0.0
        b0_z = -P.g
        #PID gains
        self.kpz = (wn_z**2 - a0_z) / b0_z
        self.kdz = (2.0*zeta_z*wn_z - a1_z) / b0_z
        self.kiz = -.01#-.2878450 #integrator gain
        self.theta_max = 10.0 *np.pi / 180.0 #maximum theta
        print('kpz = ', self.kpz)
        print('kdz = ', self.kdz)
        print('kiz = ', self.kiz)
        
        #integrator and differentiator values
        self.zdot = 0.0 #estimated derivative of z
        self.z_d1 = 0.0 #z delayed by one sample
        self.error_zdot = 0.0 #estimated derivative of error
        self.error_d1_z = 0.0 #error delayed by one sample
        self.integrator_z = 0.0 #integrator
        
    
        
    def update(self, z_r, x): 
        # extract the states
        z = x[0][0]
        theta = x[1][0]
        
        #Outer loop
        zerr = z_r - z
        #integrator and anti-windup
        self.integrator_z = self.integrator_z + P.Ts/2.0 * (zerr + self.error_d1_z)
        #differentiate z
        self.zdot = (2.0*self.sigma - P.Ts)/(2.0*self.sigma + P.Ts)*self.zdot + \
            + (2.0/(2.0*self.sigma + P.Ts))*(z - self.z_d1)
        #PID control on z
        thetaref = zerr*self.kpz + self.kiz * self.integrator_z - self.zdot*self.kdz
        thetaref_sat = self.saturate(thetaref, self.theta_max).item(0)
        # integrator anti windup
        if self.kiz != 0.0:
            self.integrator_z = self.integrator_z + P.Ts/self.kiz * (thetaref_sat - thetaref)
               
        
        #Inner loop
        thetaerr = thetaref - theta
        #differentiate theta
        self.thetadot = (2.0*self.sigma - P.Ts)/(2.0*self.sigma + P.Ts)*self.thetadot + \
            + (2.0/(2.0*self.sigma + P.Ts))*(theta - self.theta_d1)
        #? does an integrator only go on the outer loop? -> correct
        #PD control on theta 
        Ftilde_unsat = thetaerr*self.kpth - self.thetadot*self.kdth
        #saturate Ftilde
        Fe = P.m1*P.g*z/P.length + P.m2*P.g/2.0
        Fin = Ftilde_unsat + Fe
        Ftilde = self.saturate(Fin, self.force_max).item(0)
        #update delayed variables
        self.theta_d1 = theta
        self.error_d1_z = zerr
        self.z_d1 = z
        return Ftilde
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u