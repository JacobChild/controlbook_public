import numpy as np
import rodMassParam as P
import control as cnt

class ctrlObsv:
    def __init__(self):
        #tuning parameters
        tr = 1.0
        tr_obs = 0.1
        zeta = 0.707
        integrator_pole = -10.0
        zeta_obs = 0.707
        self.limit = P.tau_max
        
        #State Space Matrices
        self.A = np.array([[0.0, 1.0],
                      [-P.k1/(P.m*P.ell**2), -P.b/(P.m*P.ell**2)]])
        self.B = np.array([[0.0],
                      [1/(P.m*P.ell**2)]])
        self.C = np.array([[1.0, 0.0]])
        self.D = np.array([[0.0]])
        
        #form augmented system
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1),1)))),
                        np.hstack((-self.C, np.array([[0.0]]))) ))
        self.B1 = np.vstack((self.B, 0.0))
        #gain calculation
        wn = 2.2/tr 
        des_char_poly = np.convolve([1, 2 * zeta*wn, wn**2],
                                    [1, -integrator_pole]) #!when is the integrator pole negative vs positive?
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
        
        #? I did the below because of 3.2, but then 3.3 immediately wants a disturbance observer?
        #observer design
        wn_obs = 2.2/tr_obs 
        des_obs_char_poly = [1,2*zeta_obs*wn_obs, wn_obs**2]
        des_obs_poles = np.roots(des_obs_char_poly)
        #compute the gains if the system is observable
        if np.linalg.matrix_rank(cnt.ctrb(self.A.T, self.C.T)) != 2:
            print("The system is not observable")
        else: 
            self.L = cnt.place(self.A.T, self.C.T, des_obs_poles).T
        print('L.T: ', self.L.T)
        
        #?3.3 for disturbance observer
        #augmented matrices for observer design
        self.A2 = np.concatenate((
                            np.concatenate((self.A, self.B), axis=1),
                            np.zeros((1, 3))),
                            axis=0)
        self.B2 = np.concatenate((self.B, np.zeros((1, 1))), axis=0)
        self.C2 = np.concatenate((self.C, np.zeros((1, 1))), axis=1)
        
        #disturbance observer design
        dist_obs_pole = -1.0
        wn_obs = 2.2/tr_obs
        des_obs_char_poly = np.convolve([1, 2*zeta_obs*wn_obs, wn_obs**2],
                                        [1.0, dist_obs_pole]) #! should this pole input be negative or positive?
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
        self.thetadot = 0.0 #estimated derivative of z
        self.theta_d1 = 0.0 #z delayed by one sample
        self.integrator = 0.0
        self.error_d1 = 0.0
        self.x_hat = np.array([[0.0], #z_hat_0
                               [0.0]]) #zdot_hat_0
        self.Tau_d1 = 0.0
        self.obs_state = np.array([
            [0.0], #z_hat
            [0.0], #zdot_hat
            [0.0], # estimate of the disturbance
        ])

    def update(self, theta_r, y):
        x_hat, d_hat = self.update_observer(y)
        theta_hat = x_hat[0][0]
        error = theta_r -theta_hat
        #integrate the error
        self.integrator = self.integrator + (P.Ts/2.0)*(error + self.error_d1)
        self.error_d1 = error #update the error
        #copmute the state feedback controller
        Tau_tilde = self.K @ x_hat - self.Ki * self.integrator - d_hat
        tau = self.saturate(Tau_tilde)
        self.Tau_d1 = tau
        return tau, x_hat, d_hat

    def update_observer(self, y):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obs_state, y)
        F2 = self.observer_f(self.obs_state + P.Ts / 2 * F1, y)
        F3 = self.observer_f(self.obs_state + P.Ts / 2 * F2, y)
        F4 = self.observer_f(self.obs_state + P.Ts * F3, y)
        self.obs_state += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        x_hat = self.obs_state[0:2]
        d_hat = self.obs_state[2][0]
        return x_hat, d_hat

    def observer_f(self, x_hat, y):
        #this is called in the update observer function for RK4
        # xhat = [z_hat, zdot_hat]
        
        # xhatdot = A*(xhat-xe) + B*(u-ue) + L(y-C*xhat)
        #!is it always going to be B1 and A2 and C2 etc????
        xhat_dot = self.A2 @ x_hat\
                   + self.B1 * (self.Tau_d1)\
                   + self.L2 * (y - self.C2 @ x_hat)
        
        return xhat_dot

    def saturate(self,u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u

