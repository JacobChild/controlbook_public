import numpy as np
import control as cnt
import param as P

class Controller:
    def __init__(self):
        self.Ts = P.Ts   
        self.K = np.array([[-2.3942, -0.3995,  4.1193]])
        self.kr = 2.8021

        print('K: ', self.K)
        print('kr: ', self.kr)
        print('L^T: ', L2.T)

    def update(self, y_ref, y_m):
        x_hat, d_hat = self.update_observer(y_m)
        

        u = -self.K @ x_hat + self.kr * y_ref - d_hat
        return u.item(0), x_hat, d_hat

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
