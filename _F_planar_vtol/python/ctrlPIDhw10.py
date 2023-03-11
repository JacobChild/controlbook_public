import numpy as np
import VTOLParam as P

# set up the block beam controller class
class ctrlPD:
    def __init__(self):
        #Givens
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_h = 4.0 # rise time for altitude
        zeta_h = 0.707  # altitude damping ratio 
        M = 10.0  # Time scale separation between loops
        tr_th = 0.35 # rise time for inner loop (theta)
        zeta_th = 0.707  # inner loop damping ratio (theta)
        zeta_z = 0.707  # outer loop damping ratio (z)
        # saturation limits
        self.force_max = P.fmax  
            # maximum commanded base angle
        #---------------------------------------------------
        #                    Inner Loop
        #---------------------------------------------------
        # PD givens for Altitude (h)
        wn_h = 2.2 / tr_h
        a0_h = 0.0
        a1_h = 0.0
        b0_h = 1.0 / (P.mc + 2.0*P.mr)
        #PD gains
        self.kph = (wn_h**2 - a0_h) / b0_h
        self.kdh = (2.0*zeta_h*wn_h - a1_h) / b0_h
        self.kih = 0.05
        print('kph = ', self.kph)
        print('kdh = ', self.kdh)
        print('kih = ', self.kih)
        self.integrator_h = 0.0 #integrator
        self.h_d1 = 0.0 #h delayed by one sample
        
        # PD givens for inner loop (theta)
        wn_th = 2.2 / tr_th
        a0_th = 0.0
        a1_th = 0.0
        b0_th = 1.0 / (P.Jc + 2.0*P.mr*P.d**2)
        self.kpth = (wn_th**2 - a0_th) / b0_th
        self.kdth = (2.0*zeta_th*wn_th - a1_th) / b0_th
        print('kpth = ', self.kpth)
        print('kdth = ', self.kdth)
        
        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # PD design for outer loop (z/ roll angle phi)
        tr_z = M * tr_th  # rise time for outer loop
        wn_z =2.2 / tr_z
        a0_z = 0.0
        a1_z = P.mu / (P.mc + 2.0*P.mr)
        Fe_z = (P.mc + 2.0*P.mr)*P.g
        b0_z = -Fe_z / (P.mc + 2.0*P.mr)
        #PD gains
        self.kpz = (wn_z**2 - a0_z) / b0_z
        self.kdz = (2.0*zeta_z*wn_z - a1_z) / b0_z
        self.kiz = -.005 # integral gain
        print('kpz = ', self.kpz)
        print('kdz = ', self.kdz)
        print('kiz = ', self.kiz)
        
        #integrator and differentiator values
        self.zdot = 0.0
        self.z_d1 = 0.0
        self.error_zdot = 0.0
        self.error_z_d1 = 0.0
        self.integrator_z = 0.0
        self.sigma = 0.05
        self.thetadot = 0.0
        self.theta_d1 = 0.0
        self.error_h = 0.0
        self.hdot = 0.0
        self.error_z = 0.0
        
        
        
    def update(self, PosRef, x): 
        # extract the states -> [z, h, theta, zdot, hdot, thetadot]
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]
        
        
        #extract the reference position -> [z_r, h_r, theta_r]
        z_r = PosRef[0][0]
        h_r = PosRef[1][0]        
        
        #Longitudinal control (Altitude)
        herr = h_r - h
        #integrator and anti-windup
        if np.abs(self.hdot) < 0.01:
            self.integrator_h = self.integrator_h + (P.Ts/2.0)*(herr + self.error_h)
        else:
            self.integrator_h = 0.0 
        
        #differentiat h
        self.hdot = (2.0*self.sigma - P.Ts)/(2.0*self.sigma + P.Ts)*self.hdot + 2.0/(2.0*self.sigma + P.Ts)*(h - self.h_d1)
        
        #PID control on h
        FAltRef = self.kph*herr + self.kih*self.integrator_h - self.kdh*self.hdot
        Fe = (P.mc + 2.0*P.mr)*P.g #equilibrium force
        F = Fe + FAltRef
        F0 = self.saturate(F, self.force_max) 
        
        #Lateral control, inner and outer loop
        #Outer loop (z)
        zerr = z_r - z
        #integrator and anti-windup
        if np.abs(self.zdot) < 0.01:
            self.integrator_z = self.integrator_z + (P.Ts/2.0)*(zerr + self.error_z)
        else:
            self.integrator_z = 0.0
        
        #differentiate z
        self.zdot = (2.0*self.sigma - P.Ts)/(2.0*self.sigma + P.Ts)*self.zdot + 2.0/(2.0*self.sigma + P.Ts)*(z - self.z_d1)
        #PID control on z
        thetaref = self.kpz*zerr + self.kiz*self.integrator_z - self.kdz*self.zdot
        
        #Inner loop (theta)
        thetaerr = thetaref - theta
        #differentiate theta
        self.thetadot = (2.0*self.sigma - P.Ts)/(2.0*self.sigma + P.Ts)*self.thetadot + 2.0/(2.0*self.sigma + P.Ts)*(theta - self.theta_d1)
        #PD control on theta        
        tautilde = self.kpth*thetaerr - self.kdth*self.thetadot
        
        #Solve for fr and fl
        fl = (F0/2.0) + (tautilde/(2.0*P.d))
        fl0 = np.clip(fl, -P.fmax/2.0, P.fmax/2.0)
        #fl0 = saturate(fl, P.fmax)
        fr = tautilde/P.d + fl  
        fr0 = self.saturate(fr, P.fmax/2.0)
        #fr0 = saturate(fr, P.fmax)
        fout = np.array([[fr0], [fl0], [F0], [tautilde]])
        
        if fr > 10 or fl > 10:
            print("fr or fl was clipped!")
        
        #update the delayed variables
        self.error_z_d1 = zerr
        self.h_d1 = h
        self.theta_d1 = theta
        self.error_theta_d1 = thetaerr
        self.z_d1 = z
        self.error_z = zerr
        self.error_h = herr


        return fout
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u