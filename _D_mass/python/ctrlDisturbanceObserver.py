import numpy as np
import control as cnt
import massParam as P

class ctrlDisturbanceObserver:
    
    def __init__(self):
        #  tuning parameters
        tr = .5
        zeta = .707
        integrator_pole = -1.0 #? how do I pick the pole for integrator
        tr_obs = tr/10  # rise time frequency for observer
        zeta_obs = 0.707  # damping ratio for observer
        dist_obs_pole = 5.5 #-10.0  #? how do I pick the pole for disturbance observer
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.A = np.array([[0.0, 1.0],
                      [-P.k/P.m, -P.b/P.m]])
        self.B = np.array([[0.0],
                      [1/P.m]])
        self.C = np.array([[1.0, 0.0]])
        self.D = np.array([[0.0]])
        #form augmented system
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1),1)))),
                        np.hstack((-self.C, np.array([[0.0]]))) ))
        self.B1 = np.vstack((self.B, 0.0))
        # gain calculation
        wn = 2.2 / tr  # natural frequency
        des_char_poly = np.convolve([1, 2 * zeta*wn, wn**2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, self.B1)) != 3:
            print("The system is not controllable")
        else:
            self.K1 = (cnt.place(A1, self.B1, des_poles))
            self.K = self.K1[0][0:2]
            self.Ki = self.K1[0][2]
        print('K: ', self.K)
        print('kr: ', self.Ki)
        print(des_poles)
        
        #augmented matrices for observer design
        self.A2 = np.concatenate((
                            np.concatenate((self.A, self.B), axis=1),
                            np.zeros((1, 3))),
                            axis=0)
        self.B2 = np.concatenate((self.B, np.zeros((1, 1))), axis=0)
        self.C2 = np.concatenate((self.C, np.zeros((1, 1))), axis=1)
        
        # observer design
        wn_obs = 2.2/tr_obs
        des_obs_char_poly = np.convolve([1, 2*zeta_obs*wn_obs, wn_obs**2],
                                        [1.0, dist_obs_pole])
        des_obs_poles = np.roots(des_obs_char_poly)
        #compute the gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A2.T, self.C2.T)) != 3:
            print("The system is not observable")
        else: 
            self.L2 = cnt.place(self.A2.T, self.C2.T, des_obs_poles).T
        print('L2.T: ', self.L2.T)
        
        # dirty derivative setup
        self.sigma = 0.05 #dirty derivative gain
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.zdot = 0.0 #estimated derivative of z
        self.z_d1 = 0.0 #z delayed by one sample
        self.integrator = 0.0
        self.error_d1 = 0.0
        self.x_hat = np.array([[0.0], #z_hat_0
                               [0.0]]) #zdot_hat_0
        self.F_d1 = 0.0
        self.obs_state = np.array([
            [0.0], #z_hat
            [0.0], #zdot_hat
            [0.0], # estimate of the disturbance
        ])
        
    def update(self, z_r, zCurrent):
        # update the observer and extract z_hat
        x_hat, d_hat = self.update_observer(zCurrent)
        z_hat = x_hat[0][0]
        # extract the states
        z = zCurrent[0][0]
        error = z_r - z_hat
        #integrate the error
        self.integrator = self.integrator + (P.Ts/2.0)*(error + self.error_d1)
        self.error_d1 = error #update the error
        # compute the state feedback controller
        F_tilde = -self.K @ x_hat -self.Ki * self.integrator - d_hat
        Fout = self.saturate(F_tilde.item(0), P.F_max)
        #print(Fout, "Fout")
        self.F_d1 = Fout 
        return Fout, x_hat, d_hat
    
    def update_observer(self, z_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obs_state, z_m)
        F2 = self.observer_f(self.obs_state + P.Ts / 2 * F1, z_m)
        F3 = self.observer_f(self.obs_state + P.Ts / 2 * F2, z_m)
        F4 = self.observer_f(self.obs_state + P.Ts * F3, z_m)
        self.obs_state += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        x_hat = self.obs_state[0:2]
        d_hat = self.obs_state[2][0]
        return x_hat, d_hat

    def observer_f(self, x_hat, z_m):
        # xhat = [z_hat, zdot_hat]
        
        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        xhat_dot = self.A2 @ x_hat\
                   + self.B1 * (self.F_d1)\
                   + self.L2 * (z_m - self.C2 @ x_hat)
        return xhat_dot
        
    
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u