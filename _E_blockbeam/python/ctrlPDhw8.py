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
        tr_th = 0.50 # rise time for inner loop
        zeta_th = 0.707  # inner loop damping ratio 
        M = 10.0  # Time scale separation between loops
        zeta_z = 0.707  # outer loop damping ratio
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
        #PD gains
        self.kpth = (wn_th**2 - a0_th) / b0_th
        self.kdth = (2.0*zeta_th*wn_th - a1_th) / b0_th
        print('kpth = ', self.kpth)
        print('kdth = ', self.kdth)
        
        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # PD design for outer loop
        tr_z = M * tr_th  # rise time for outer loop
        wn_z =2.2 / tr_z
        a0_z = 0.0
        a1_z = 0.0
        z0_err_z = 0.0
        b0_z = -P.g
        #PD gains
        self.kpz = (wn_z**2 - a0_z) / b0_z
        self.kdz = (2.0*zeta_z*wn_z - a1_z) / b0_z
        print('kpz = ', self.kpz)
        print('kdz = ', self.kdz)
        
    def update(self, z_r, x): 
        # extract the states
        z = x[0][0]
        theta = x[1][0]
        zdot = x[2][0]
        thetadot = x[3][0]
        
        #Outer loop
        thetaref = (z_r-z)*self.kpz - zdot*self.kdz
        #Inner loop
        Ftilde = self.kpth*(thetaref-theta) - self.kdth*thetadot
        Fe = P.m1*P.g*z/P.length + P.m2*P.g/2.0
        Fnew = Ftilde + Fe
        return saturate(Fnew, P.Fmax)
    
def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u