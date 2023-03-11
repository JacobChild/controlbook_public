import numpy as np
import control as cnt
import VTOLParam as P

#set up block beam controller class
class ctrlStateFeedback:
    def __init__(self):
        #tuning parameters
        tr_h = 3.0 # rise time for altitude
        zeta_h = 0.707  # altitude damping ratio 
        M = 10.0  # Time scale separation between loops
        tr_th = 0.45 # rise time for inner loop (theta)
        tr_z = M * tr_th  # rise time for outer loop
        zeta_th = 0.707  # inner loop damping ratio (theta)
        zeta_z = 0.707  # outer loop damping ratio (z)
        wn_h = 2.2 / tr_h
        wn_th = 2.2 / tr_th
        wn_z = 2.2 / tr_z
        # State Space Equations
        # xdot = Ax + Bu
        # y = Cx + Du
        
        #Longitudinal State Space
        Alon = np.array([[0.0, 1.0],
                         [0.0, 0.0]])
        Blon = np.array([[ 0.0],
                         [1.0/ (P.mc + 2.0*P.mr)]])
        Clon = np.array([[1.0, 0.0]])
        Dlon = np.array([[0.0]])
        des_char_poly_lon = [1.0, 2.0*zeta_h*wn_h, wn_h**2]
        des_poles_lon = np.roots(des_char_poly_lon)
        print(des_poles_lon, "longitudinal poles")
        
        #compute the longitudinal gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(Alon, Blon)) != 2:
            print("The longitudinal system is not controllable")
        else: 
            self.Klon = (cnt.acker(Alon, Blon, des_poles_lon))
            self.krlon = -1.0 / (Clon @ np.linalg.inv(Alon - Blon @ self.Klon) @ Blon)
        print('Klon: ', self.Klon)
        print('krlon: ', self.krlon)
        
        #Lateral State Space
        Fe = (P.mc + 2.0*P.mr)*P.g
        Alatterm1 = -Fe / (P.mc + 2.0*P.mr)
        Alatterm2 = -P.mu / (P.mc + 2.0*P.mr)
        Blatterm1 = 1.0 / (P.Jc + 2.0*P.mr*P.d**2)
        Alat = np.array([[0.0, 0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0, 1.0],
                         [0.0, Alatterm1, Alatterm2, 0.0],
                         [0.0, 0.0, 0.0, 0.0]])
        Blat = np.array([[0.0],
                         [0.0],
                         [0.0],
                         [Blatterm1]])
        Clat = np.array([[1.0, 0.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0, 0.0]])
        Dlat = np.array([[0.0],
                         [0.0]])
        #lateral gain calculation
        des_char_poly_lat = np.convolve([1, 2.0*zeta_th*wn_th, wn_th**2],
                                        [1, 2.0*zeta_z*wn_z, wn_z**2])
        des_poles_lat = np.roots(des_char_poly_lat)

        #compute the lateral gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(Alat, Blat)) != 4:
            print("The lateral system is not controllable")
        else: 
            self.Klat = (cnt.acker(Alat, Blat, des_poles_lat))
            Crlat = np.array([[1.0, 0.0, 0.0, 0.0]])
            self.krlat = -1.0 / (Crlat @ np.linalg.inv(Alat - Blat @ self.Klat) @ Blat)
            #?self.krlat = -1.0 / (Clat @ np.linalg.inv(Alat - Blat @ self.Klat) @ Blat)
            print('Klat: ', self.Klat)
            print('krlat: ', self.krlat)
            
            
    def update(self, PosRef, x):
        # extract the states
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]
        zdot = x[3][0]
        hdot = x[4][0]
        thetadot = x[5][0]
        LonStates = np.array([[h], [hdot]])
        LatStates = np.array([[z], [theta], [zdot], [thetadot]])
        
        # extract reference position
        z_r = PosRef[0][0]
        h_r = PosRef[1][0]
        
        #Longitudinal Control
        Fe = (P.mc + 2.0*P.mr)*P.g #equilibrium force
        Ftilde = -self.Klon @ LonStates + self.krlon * h_r 
        Flon = Ftilde + Fe
        F0 = self.saturate(Flon, P.fmax)
        
    
        #Lateral Control
        tautilde = -self.Klat @ LatStates  + self.krlat * z_r
        
        #Solve for fr and fl
        fl = (F0.item(0)/2.0) + (tautilde.item(0)/(2.0*P.d))
        fl0 = np.clip(fl, -P.fmax/2.0, P.fmax/2.0)
        #fl0 = saturate(fl, P.fmax)
        fr = tautilde.item(0)/P.d + fl  
        fr0 = self.saturate(fr, P.fmax/2.0)
        #fr0 = saturate(fr, P.fmax)
        fout = np.array([[fr0], [fl0], [F0.item(0)], [tautilde.item(0)]])
        
        return fout
        
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u