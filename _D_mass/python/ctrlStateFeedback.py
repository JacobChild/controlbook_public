import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedback:
    
    def __init__(self):
        #  tuning parameters
        tr = 1.0
        zeta = .707
        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([[0.0, 1.0],
                      [-P.k/P.m, -P.b/P.m]])
        B = np.array([[0.0],
                      [1/P.m]])
        C = np.array([[1.0, 0.0]])
        D = np.array([[0.0]])
        # gain calculation
        wn = 2.2 / tr  # natural frequency
        des_char_poly = [1, 2 * zeta*wn, wn**2]
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 2:
            print("The system is not controllable")
        else:
            self.K = (cnt.acker(A, B, des_poles))
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print(des_poles)
        # dirty derivative setup
        self.sigma = 0.05 #dirty derivative gain
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.zdot = 0.0 #estimated derivative of z
        self.z_d1 = 0.0 #z delayed by one sample
        
    def update(self, z_r, zCurrent):
        z = zCurrent[0][0]
        zdotstate = zCurrent[1][0]
        # differentiate z to get zdot using the dirty derivative
        self.zdot = self.beta*self.zdot + (1.0-self.beta)*((z - self.z_d1)/P.Ts)
        x = np.array([[z], [self.zdot]])
        F_tilde = -self.K @ x + self.kr * z_r
        Fout = self.saturate(F_tilde, P.F_max)
        #print(Fout, "Fout")
        self.z_d1 = z 
        return Fout.item(0)
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u