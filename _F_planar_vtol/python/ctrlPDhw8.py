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
        print('kph = ', self.kph)
        print('kdh = ', self.kdh)
        
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
        print('kpz = ', self.kpz)
        print('kdz = ', self.kdz)
        
    def update(self, PosRef, x): 
        # extract the states -> [z, h, theta, zdot, hdot, thetadot]
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]
        zdot = x[3][0]
        hdot = x[4][0]
        thetadot = x[5][0]
        
        #extract the reference position -> [z_r, h_r, theta_r]
        z_r = PosRef[0][0]
        h_r = PosRef[1][0]        
        
        #Longitudinal control (Altitude)
        FAltRef = self.kph*(h_r-h) - self.kdh*hdot
        Fe = (P.mc + 2.0*P.mr)*P.g #! when do I use Fe vs not?
        F = Fe + FAltRef
        F0 = np.clip(F, 0.0, P.fmax) 
        
        #Lateral control, inner and outer loop
        #Outer loop (z)
        thetaref = self.kpz*(z_r-z) - self.kdz*zdot
        
        #Inner loop (theta)
        tautilde = self.kpth*(thetaref-theta) - self.kdth*thetadot
        
        #Solve for fr and fl
        fl = (F0/2.0) + (tautilde/(2.0*P.d))
        fl0 = np.clip(fl, -P.fmax/2.0, P.fmax/2.0)
        #fl0 = saturate(fl, P.fmax)
        fr = tautilde/P.d + fl  
        fr0 = np.clip(fr, -P.fmax/2.0, P.fmax/2.0)
        #fr0 = saturate(fr, P.fmax)
        fout = np.array([[fr0], [fl0], [F0], [tautilde]])
        
        if fr > 10 or fl > 10:
            print("fr or fl was clipped!")
        
        
        return fout
    
def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u