import numpy as np
import blockbeamParam as P
import control as cnt

# set up the block beam controller class
class ctrlDisturbanceObserver:
    def __init__(self):
        #-----
        # State Feedback Control Design
        #----
        # tuning parameters
        tr_z = 1.5 #1.2 # rise time for outer loop
        zeta_z = 0.707  # inner loop damping ratio
        M = 10.0  # Time scale separation between loops
        tr_th = tr_z/M #0.5  # rise time for inner loop
        zeta_th = 0.707  # outer loop damping ratio
        tr_z_obs = tr_z/5.0 # rise time for position
        tr_theta_obs = tr_th / 5.0  # rise time for angle
        self.sigma = 0.05  # cutoff freq for dirty derivative
        # saturation limits
        self.force_max = P.Fmax  
        wn_th = 2.2 / tr_th
        z0_err_th = P.length/2.0
        b0_th = P.length/(P.m2*P.length**2 / 3.0 + P.m1 * z0_err_th**2)
        wn_z =2.2 / tr_z
        integrator_pole = -4.0
        dist_obsv_pole = -7.25 #? how do you pick the poles?
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.Ze = P.length/2.0
        Aterm1 = -P.m1*P.g/((P.m2*P.length**2)/3.0 + P.m1*self.Ze**2)
        self.A = np.array([[0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [0.0, -P.g, 0.0, 0.0],
                      [Aterm1, 0.0, 0.0, 0.0]])
        Bterm1 = P.length/(P.m2*P.length**2/3.0 + P.m1*self.Ze**2)
        self.B = np.array([[0.0], [0.0], [0.0], [Bterm1]])
        self.C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        D = np.array([[0.0]])
        # form augmented system
        self.Cr = np.array([[1, 0]]) @ self.C
        A1 = np.concatenate((
                np.concatenate((self.A, np.zeros((4, 1))), axis=1),
                np.concatenate((-self.Cr, np.matrix([[0.0]])), axis=1)),
                axis=0)
        self.B1 = np.concatenate((self.B, np.matrix([[0.0]])), axis=0)
        # control gain calculation
        des_char_poly = np.convolve(
                np.convolve([1, 2 * zeta_z * wn_z, wn_z**2],
                            [1, 2 * zeta_th * wn_th, wn_th**2]),
                np.poly([integrator_pole]))
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, self.B1)) != 5:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, self.B1, des_poles)
            self.K = K1[0][0:4]
            self.ki = K1[0][4]
        print('K: ', self.K)
        print('ki: ', self.ki)
        #compute observer gains
        #Aurmented matrices
        self.A2 = np.concatenate((
            np.concatenate((self.A, self.B), axis=1),
            np.zeros((1, 5))), axis=0)
        self.C2 = np.concatenate((self.C, np.zeros((2, 1))), axis=1)
        
        wn_z_obs = 2.2/tr_z_obs
        wn_th_obs = 2.2/tr_theta_obs
        des_obsv_char_poly = np.convolve(
            np.convolve([1, 2*zeta_z*wn_z_obs, wn_z_obs**2],
            [1, 2*zeta_th*wn_th_obs, wn_th_obs**2]),
            [1.0, -dist_obsv_pole])
        des_obs_poles = np.roots(des_obsv_char_poly)
        #compute the observer gain if system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A2.T, self.C2.T)) != 5:
            print("The system is not observable")
        else:
            self.L2 = cnt.place(self.A2.T, self.C2.T, des_obs_poles).T 
        print('L2.T: ', self.L2.T)               
        # variables to implement integrator
        self.integrator_z = 0.0 #integrator
        self.theta_max = 10.0 *np.pi / 180.0 #maximum theta
        self.error_d1_z = 0.0 #error delayed by one sample
        # estimated state variables
        self.observer_state = np.array([
            [0.0],
            [0.0],
            [0.0],
            [0.0],
            [0.0]])
        self.F_d1 = 0.0 # delayed force
        
        
    
        
    def update(self, z_r, x): 
        # update the observer and extract z_hat
        #print("x: ", x)
        x_hat, d_hat = self.update_observer(x)
        z_hat = self.Cr @ x_hat
        zerror = z_r - z_hat
        self.integrator_z = self.integrator_z + (P.Ts/2.0)*(zerror + self.error_d1_z)
        self.error_d1_z = zerror
        #construct the linearized state
        xe = np.array([[P.z0], [0.0], [0.0], [0.0]])
        x_tilde = x_hat - xe
        Ftilde = -self.K @ x_tilde - self.ki*self.integrator_z 
        Fe = P.m1*P.g*self.Ze/P.length + P.m2*P.g/2.0 
        Fin = Ftilde[0][0] + Fe - d_hat
        Fout = self.saturate(Fin, P.Fmax)
        # integrator anti-windup
        if self.ki != 0.0:
            self.integrator_z = self.integrator_z + P.Ts/self.ki * (Fout - Fin)
        self.F_d1 = Fout.item(0)
        return Fout.item(0), x_hat, d_hat
    
    def update_observer(self, z_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.observer_state, z_m)
        F2 = self.observer_f(self.observer_state + P.Ts / 2 * F1, z_m)
        F3 = self.observer_f(self.observer_state + P.Ts / 2 * F2, z_m)
        F4 = self.observer_f(self.observer_state + P.Ts * F3, z_m)
        self.observer_state += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        x_hat = self.observer_state[0:4]
        d_hat = self.observer_state[4][0]
        return x_hat, d_hat
    
    def observer_f(self, x_hat, z_m):
        # xhat = [z_hat, zdot_hat] #? should xhat be z, theta, zdot, thetadot? predicitons?
        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        xe = np.array([[self.Ze], [0.0], [0.0], [0.0], [0.0]])
        Fe = P.m1*P.g*self.Ze/P.length + P.m2*P.g/2.0
        xhat_dot = self.A2 @ (x_hat-xe)\
                   + self.B1 * (self.F_d1- Fe)\
                   + self.L2 @ (z_m - self.C2 @ x_hat)
        return xhat_dot
    
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u