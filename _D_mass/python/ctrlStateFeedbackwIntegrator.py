import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedback:
    
    def __init__(self):
        #  tuning parameters
        tr = 1.0
        zeta = .707
        integrator_pole = -2.0
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([[0.0, 1.0],
                      [-P.k/P.m, -P.b/P.m]])
        B = np.array([[0.0],
                      [1/P.m]])
        C = np.array([[1.0, 0.0]])
        D = np.array([[0.0]])
        #form augmented system
        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A, 1),1)))),
                        np.hstack((-C, np.array([[0.0]]))) ))
        B1 = np.vstack((B, 0.0))
        # gain calculation
        wn = 2.2 / tr  # natural frequency
        des_char_poly = np.convolve([1, 2 * zeta*wn, wn**2],
                                    [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 3:
            print("The system is not controllable")
        else:
            self.K1 = (cnt.place(A1, B1, des_poles))
            self.K = self.K1[0][0:2]
            self.Ki = self.K1[0][2]
        print('K: ', self.K)
        print('kr: ', self.Ki)
        print(des_poles)
        # dirty derivative setup
        self.sigma = 0.05 #dirty derivative gain
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.zdot = 0.0 #estimated derivative of z
        self.z_d1 = 0.0 #z delayed by one sample
        self.integrator = 0.0
        self.error_d1 = 0.0
        
    def update(self, z_r, zCurrent):
        z = zCurrent[0][0]
        zdotstate = zCurrent[1][0]
        error = z_r - z
        #integrate the error
        self.integrator = self.integrator + (P.Ts/2.0)*(error + self.error_d1)
        self.error_d1 = error #update the error
        # differentiate z to get zdot using the dirty derivative
        self.zdot = self.beta*self.zdot + (1.0-self.beta)*((z - self.z_d1)/P.Ts)
        x = np.array([[z], [self.zdot]])
        F_tilde = -self.K @ x -self.Ki * self.integrator
        Fout = self.saturate(F_tilde, P.F_max)
        #print(Fout, "Fout")
        self.z_d1 = z 
        return Fout.item(0)
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u