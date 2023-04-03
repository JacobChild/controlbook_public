import numpy as np
import control as cnt
import param as P

#! I know this isn't good, but I don't have time


class Controller:
    def __init__(self):
        tr1 = 1.0
        tr2 = tr1
        tr3 = tr1 / 10.0
        zeta = 0.707
        wn1 = 2.2 / tr1
        wn2 = 2.2 / tr2
        wn3 = 2.2 / tr3
        self.integrator1poles = [-1.0]
        self.integrator2poles = [-1.0]
        self.A = np.array([[-1.0, 2.0, -3.0],
                      [4.0, -5.0, 6.0],
                      [0.0, 1.0, 0.0]])
        self.B = np.array([[7.0], [-8.0], [9.0]])
        self.C = np.array([1.0, -1.0, 1.0])
        des_char_poly = np.convolve(np.convolve([1.0, 2.0*zeta*wn1, wn1**2],
                                    [1.0, 2.0*zeta*wn2, wn2**2]),
                                    [1.0, -self.integrator1poles[0]])
        des_poles_1 = np.roots(des_char_poly)
        # form augmented system
        #form augmented system
        A1 = np.vstack((np.hstack((self.A, np.zeros((np.size(self.A, 1),1)))),
                        np.hstack((-self.C, np.array([[0.0]]))) ))
        B1 = np.vstack((self.B, 0.0))
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != 4:
            print("The system is not controllable")
        else:
            self.K1 = (cnt.place(A1, B1, des_poles_1))
            self.K = self.K1[0][0:3]
            self.Ki = self.K1[0][3]
        print('K: ', self.K)
        print('ki: ', self.Ki)
        
        self.integrator1 = 0.0
        self.error_d1_1 = 0.0

    def update(self, y_ref, x):
        y = P.C @ x
        #create the states
        
        error = y_ref - y[0][0]
        self.integrator1 = self.integrator1 + (P.Ts/2.0)*(error + self.error_d1_1)
        self.error_d1_1 = error
        u = -self.K @ y - self.Ki * self.integrator1
        return u


