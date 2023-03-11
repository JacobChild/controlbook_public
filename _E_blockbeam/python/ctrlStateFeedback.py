import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedback:
    
    def __init__(self):
        tr_z = 1.0
        zeta_z = 0.707
        M = 10.0 #time scale separation
        tr_th = tr_z / M
        zeta_th = 0.707
        self.sigma = 0.05
        wn_z = 2.2 / tr_z
        wn_th = 2.2 / tr_th
        Ze = P.length/2.0
        Aterm1 = -P.g * P.m1 / (P.m2*P.length**2/3.0 + P.m1 * Ze**2 )
        Bterm1 = P.length / (P.m2*P.length**2/3.0 + P.m1 * Ze**2)
        A = np.array([[0.0, 0.0, 1.0, 0.0],
                      [0.0, 0.0, 0.0, 1.0],
                      [0.0, -P.g, 0.0, 0.0],
                      [Aterm1, 0.0, 0.0, 0.0]])
        B = np.array([[0.0],
                      [0.0],
                      [0.0],
                      [Bterm1]])
        C = np.array([[1.0, 0.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0, 0.0]])
        D = np.array([[0.0]])
                     
        des_char_poly = np.convolve([1, 2*zeta_z*wn_z, wn_z**2],
                                    [1, 2*zeta_th*wn_th, wn_th**2])
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != 4:
            print("The system is not controllable")
        else:
            self.K = cnt.acker(A, B, des_poles)
            Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.kr = -1.0 / (Cr @ np.linalg.inv(A - B @ self.K) @ B)
        # print gains to terminal
        print('K: ', self.K)
        print('kr: ', self.kr)
        
        #? no dirty derivative setup?
        
    def update(self, z_r, x):
        z = x[0][0]
        #print(self.K, "self.K")
        #print(x, "x")
        Ftilde = -self.K @ x + self.kr * z_r
        Fe = P.m1*P.g*z/P.length + P.m2*P.g/2.0
        Fin = Ftilde + Fe
        Fout = self.saturate(Fin, P.Fmax)
        return Fout.item(0)
    
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u