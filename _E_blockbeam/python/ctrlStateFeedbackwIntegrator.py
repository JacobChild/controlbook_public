import numpy as np
import blockbeamParam as P
import control as cnt

# set up the block beam controller class
class ctrlStateFeedbackwIntegrator:
    def __init__(self):
        #-----
        # State Feedback Control Design
        #----
        # tuning parameters
        tr_z = 1.5 #1.2 # rise time for outer loop
        zeta_th = 0.707  # inner loop damping ratio 
        M = 10.0  # Time scale separation between loops
        tr_th = tr_z/M #0.5  # rise time for inner loop
        zeta_z = 0.707  # outer loop damping ratio
        self.sigma = 0.05  # cutoff freq for dirty derivative
        # saturation limits
        self.force_max = P.Fmax  
        wn_th = 2.2 / tr_th
        z0_err_th = P.length/2.0
        b0_th = P.length/(P.m2*P.length**2 / 3.0 + P.m1 * z0_err_th**2)
        wn_z =2.2 / tr_z
        integrator_pole = -1.0
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        Ze = P.length/2.0
        Aterm1 = -P.m1*P.g/((P.m2*P.length**2)/3.0 + P.m1*Ze**2)
        A = np.array([[0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [0.0, -P.g, 0.0, 0.0],
                      [Aterm1, 0.0, 0.0, 0.0]])
        Bterm1 = P.length/(P.m2*P.length**2/3.0 + P.m1*Ze**2)
        B = np.array([[0.0], [0.0], [0.0], [Bterm1]])
        C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        D = np.array([[0.0]])
        # form augmented system
        self.Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
        A1 = np.vstack((np.hstack((A, np.zeros((4,1)))),
                        np.hstack((-self.Cr, np.zeros((1,1)) )) ))
        B1 = np.vstack((B, np.zeros((1,1)) ))
        # gain calculation
        des_char_poly = np.convolve(
            np.convolve([1, 2*zeta_z*wn_z, wn_z**2],
                        [1, 2*zeta_th*wn_th, wn_th**2]),
            [1.0, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        print('K: ', self.K)
        print('ki: ', self.ki)
        # variables to implement integrator
    
        #self.thetadot = 0.0 #estimated derivative of theta
        #self.theta_d1 = 0.0 #theta delayed by one sample
        #self.error_thetadot = 0.0 #estimated derivative of error
        self.integrator_z = 0.0 #integrator
        self.theta_max = 10.0 *np.pi / 180.0 #maximum theta
        #self.zdot = 0.0 #estimated derivative of z
        #self.z_d1 = 0.0 #z delayed by one sample
        #self.error_zdot = 0.0 #estimated derivative of error
        self.error_d1_z = 0.0 #error delayed by one sample
        
    
        
    def update(self, z_r, x): 
        z = self.Cr @ x
        #print(self.K, "self.K")
        #print(x, "x")
        # create the current state
        x_tilde = x - np.array([[P.z0], [0.0], [0.0], [0.0]])
        zerror = z_r - z
        self.integrator_z = self.integrator_z + (P.Ts/2.0)*(zerror + self.error_d1_z)
        self.error_d1_z = zerror
        Ftilde = -self.K @ x_tilde - self.ki*self.integrator_z 
        Fe = P.m1*P.g*P.z0/P.length + P.m2*P.g/2.0 
        Fin = Ftilde + Fe
        Fout = self.saturate(Fin, P.Fmax)
        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator_z = self.integrator_z + P.Ts/self.ki * (Fout - Fin)
            
        return Fout.item(0)
        
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u