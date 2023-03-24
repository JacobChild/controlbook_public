import numpy as np
import control as cnt
import VTOLParam as P

#set up block beam controller class
class ctrlStateFeedback:
    def __init__(self):
        #tuning parameters
        tr_h = 1.0 # rise time for altitude
        zeta_h = 0.707  # altitude damping ratio 
        M = 10.0  # Time scale separation between loops
        tr_th = 0.2 # rise time for inner loop (theta)
        tr_z = M * tr_th  # rise time for outer loop
        zeta_th = 0.707  # inner loop damping ratio (theta)
        zeta_z = 0.707  # outer loop damping ratio (z)
        wn_h = 2.2 / tr_h
        wn_th = 2.2 / tr_th
        wn_z = 2.2 / tr_z
        self.integrator_h_poles = [-1.0]  # integrator on h -> helps with poles?
        self.integrator_z_poles = [-1.0]  # integrator on z
        # State Space Equations
        # xdot = Ax + Bu
        # y = Cx + Du
        
        #Longitudinal State Space, hard coded rather than np.hstack and np.vstack
        Alon = np.array([[0.0, 1.0, 0.0],
                         [0.0, 0.0, 0.0],
                         [-1.0, 0.0, 0.0]])
        Blon = np.array([[ 0.0],
                         [1.0/ (P.mc + 2.0*P.mr)],
                         [0.0]])
        Clon = np.array([[1.0, 0.0, 0.0]])
        Dlon = np.array([[0.0]])
        des_char_poly_lon = np.convolve([1.0, 2.0*zeta_h*wn_h, wn_h**2],
                                        np.poly(self.integrator_h_poles))
        des_poles_lon = np.roots(des_char_poly_lon)
        print(des_poles_lon, "longitudinal poles")
        
        #compute the longitudinal gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(Alon, Blon)) != 3:
            print("The longitudinal system is not controllable")
        else: 
            self.K1lon = (cnt.place(Alon, Blon, des_poles_lon))
            self.Klon = np.array([self.K1lon.item(0), self.K1lon.item(1)]) 
            self.ki_lon = self.K1lon.item(2)
        print('Klon: ', self.Klon)
        print('ki_lon: ', self.ki_lon)
        
        #Lateral State Space
        Fe = (P.mc + 2.0*P.mr)*P.g
        Alatterm1 = -Fe / (P.mc + 2.0*P.mr)
        Alatterm2 = -P.mu / (P.mc + 2.0*P.mr)
        Blatterm1 = 1.0 / (P.Jc + 2.0*P.mr*P.d**2)
        Alat = np.array([[0.0, 0.0, 1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 1.0, 0.0],
                         [0.0, Alatterm1, Alatterm2, 0.0, 0.0],
                         [0.0, 0.0, 0.0, 0.0, 0.0],
                         [-1.0, 0.0, 0.0, 0.0, 0.0]])
        Blat = np.array([[0.0],
                         [0.0],
                         [0.0],
                         [Blatterm1],
                         [0.0]])
        Clat = np.array([[1.0, 0.0, 0.0, 0.0, 0.0],
                         [0.0, 1.0, 0.0, 0.0, 0.0]])
        Dlat = np.array([[0.0],
                         [0.0]])
        Cr = np.array([[1.0, 0.0, 0.0, 0.0]])
        #A1 = np.vstack((np.hstack((Alat, np.zeros((np.size(Alat,1),1)))),
                        #np.hstack((-Cr, np.array([[0.0]]))) ))
        #B1 = np.vstack( (Blat, 0.0) )
        #lateral gain calculation
        des_char_poly_lat = np.convolve(
            np.convolve([1, 2.0*zeta_z*wn_z, wn_z**2],
                        [1, 2.0*zeta_th*wn_th, wn_th**2]),
            np.poly(self.integrator_z_poles))
        des_poles_lat = np.roots(des_char_poly_lat)

        #compute the lateral gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(Alat, Blat)) != 5:
            print("The lateral system is not controllable")
        else: 
            K1 = (cnt.place(Alat, Blat, des_poles_lat))
            self.K_lat = K1[0][0:4]
            self.ki_lat = K1[0][4]
            print('K_lat: ', self.K_lat)
            print('ki_lat: ', self.ki_lat)
        # variables to implement integrator
        self.integrator_z = 0.0
        self.error_z_d1 = 0.0
        self.z_d1 = 0.0
        self.integrator_h = 0.0
        self.error_h_d1 = 0.0
        self.h_d1 = 0.0
        self.limit = P.fmax
        self.F_limit = P.fmax * 2.0
        self.tau_limit = P.fmax * P.d * 2.0
        self.sigma = 0.05
        self.beta = (2.0 * self.sigma - P.Ts) / (2.0 * self.sigma + P.Ts) #dirty derivative gain
        self.z_dot = 0.0
        self.h_dot = 0.0
        self.theta_dot = 0.0
        self.theta_d1 = 0.0
        
            
            
    def update(self, PosRef, x):
        # extract the states
        z = x[0][0]
        h = x[1][0]
        theta = x[2][0]

        # Construct the states
        LonStates = np.array([[h], [self.differentiateH(h)]])
        LatStates = np.array([[z], [theta], [self.differentiateZ(z)], [self.differentiateTheta(theta)]])
        
        # extract reference position
        z_r = PosRef[0][0]
        h_r = PosRef[1][0]
        error_z = z_r - z
        error_h = h_r - h
        
        # integrate error
        self.integrateErrorZ(error_z)
        self.integrateErrorH(error_h)
        
        # Compute the state feedback controllers
        
        #Longitudinal Control
        Fe = (P.mc + 2.0*P.mr)*P.g/np.cos(theta) #equilibrium force
        Ftilde = -self.Klon @ LonStates -self.ki_lon * self.integrator_h
        Flon = Ftilde + Fe
        Fout = self.saturate(Flon.item(0), P.fmax)
        #print("Fout:", Fout)
        self.integratorAntiWindup(Fout, Flon, self.ki_lon, self.integrator_h)
        
        #Lateral Control
        tautilde = -self.K_lat @ LatStates  - self.ki_lat * self.integrator_z
        tausat = self.saturate(tautilde.item(0), self.tau_limit)
        self.integratorAntiWindup(tausat, tautilde, self.ki_lat, self.integrator_z)
        
        #Solve for fr and fl
        #fl = (Fout.item(0)/2.0) + (tautilde.item(0)/(2.0*P.d))
        #fl0 = self.saturate(fl, P.fmax/2.0) #np.clip(fl, -P.fmax/2.0, P.fmax/2.0)
        #fl0 = saturate(fl, P.fmax)
        #fr = (Fout.item(0)/2.0) - (tautilde.item(0)/(2.0*P.d))  
        #fr0 = self.saturate(fr, P.fmax/2.0)
        #fr0 = saturate(fr, P.fmax)
        #fout1 = np.array([[fr0], [fl0], [Fout.item(0)], [tausat.item(0)]])
        fout1 = np.array([[Fout], [tausat]])
        
        return fout1
    
    # Extra functions to differentiate, integrate, anti windup, and saturate
    
    def differentiateZ(self,z):
        self.z_dot = self.beta*self.z_dot + (1-self.beta)*((z-self.z_d1)/P.Ts)
        self.z_d1 = z
        return self.z_dot
        
    def differentiateH(self,h):
        self.h_dot = self.beta*self.h_dot + (1-self.beta)*((h-self.h_d1)/P.Ts)
        self.h_d1 = h
        return self.h_dot
        
    def differentiateTheta(self, theta):
        self.theta_dot = self.beta*self.theta_dot + (1-self.beta)*((theta-self.theta_d1)/P.Ts)
        self.theta_d1 = theta
        return self.theta_dot
        
    def integrateErrorZ(self, error_z):
        self.integrator_z = self.integrator_z + (P.Ts/2.0)*(error_z + self.error_z_d1)
        self.error_z_d1 = error_z
        
    def integrateErrorH(self, error_h):
        self.integrator_h = self.integrator_h + (P.Ts/2.0)*(error_h + self.error_h_d1)
        self.error_h_d1 = error_h
        
    def integratorAntiWindup(self, u_sat, u_unsat, ki, integrator):
        if ki != 0.0:
            integrator = integrator + P.Ts/ki * (u_sat - u_unsat)   
        
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
            #print("Saturation")
        return u