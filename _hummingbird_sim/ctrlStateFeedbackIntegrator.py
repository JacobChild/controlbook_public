import numpy as np
import control as cnt
import hummingbirdParam as P
import sympy as sym


class ctrlStateFeedbackIntegrator:
    def __init__(self):
        
        
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        # tuning parameters
        wn_th = 2.4
        zeta_th = .707  
        pi_th = 4.0
        wn_psi = 1.6
        zeta_psi = .707
        wn_phi = 10.0
        zeta_phi = .707  
        pi_psi = 1.9
        
        # soft code
        #Set up to compare acker vs hardcoded
        JT = P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J2z + P.m3 * (P.ell3x**2 + P.ell3y**2)
        Fe = (P.m1*P.ell1 + P.m2*P.ell2)*P.g / P.ellT
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        print('b_theta',b_theta)
        ALon = np.array([[0.0,1.0],
                        [0.0,0.0]])
        BLon = np.array([[0.0],[b_theta]])
        CrLon = np.array([[1.0,0.0]])
        ALatTerm1 = P.ellT * Fe / (JT + P.J1z)
        print('ALatTerm1',ALatTerm1)
        ALat = np.array([[0.0,0.0,1.0,0.0],
                        [0.0,0.0,0.0,1.0],
                        [0.0,0.0,0.0,0.0],
                        [ALatTerm1,0.0,0.0,0.0]])
        BLat = np.array([[0.0],
                        [0.0],
                        [1.0/(P.J1x)],
                        [0.0]])
        CrLat = np.array([[0.0,1.0,0.0,0.0]])

        ALonAug = np.vstack((np.hstack((ALon, np.zeros((np.size(ALon,1),1)))),
                                    np.hstack((-CrLon, np.array([[0.0]]))) ))
        BLonAug = np.vstack((BLon, 0.0))
        print('ALonAug',ALonAug)
        print('BLonAug',BLonAug)
        ALatAug = np.vstack((np.hstack((ALat, np.zeros((np.size(ALat,1),1)))),
                            np.hstack((-CrLat, np.array([[0.0]]))) ))
        BLatAug = np.vstack((BLat, 0.0))
        print('ALatAug',ALatAug)
        print('BLatAug',BLatAug)
        des_char_poly_lon = np.convolve([1.0, 2.0*zeta_th*wn_th, wn_th**2],
                                        np.poly([-4.0]))
        des_poles_lon = np.roots(des_char_poly_lon)
        if np.linalg.matrix_rank(cnt.ctrb(ALonAug, BLonAug)) != 3:
            print("The longitudinal system is not controllable")
        else:
            self.K1lon = cnt.acker(ALonAug, BLonAug, des_poles_lon)
            
        #Lateral State space
        des_char_poly_lat = np.convolve(np.poly([-1.9]),
                        np.convolve([1.0, 2.0*zeta_psi*wn_psi, wn_psi**2],
                        [1.0, 2.0*zeta_phi*wn_phi, wn_phi**2]))
        des_poles_lat = np.roots(des_char_poly_lat)
        if np.linalg.matrix_rank(cnt.ctrb(ALatAug, BLatAug)) != 5:
            print("The lateral system is not controllable")
        else:
            self.K1lat = cnt.place(ALatAug, BLatAug, des_poles_lat)
            
        print("des_char_poly_lon: ", des_char_poly_lon)
        
        print("des_char_poly_lat: ", des_char_poly_lat)
        print('K1lon: ', self.K1lon)
        print('K1lat: ', self.K1lat)
        # hard code Ackerman's formula
        JT = P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J2z + P.m3 * (P.ell3x**2 + P.ell3y**2)
        Fe = (P.m1*P.ell1 + P.m2*P.ell2)*P.g / P.ellT
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        alpha1_lon = pi_th + 2*wn_th*zeta_th
        alpha2_lon = 2*pi_th*wn_th*zeta_th + wn_th**2
        alpha3_lon = pi_th*wn_th**2 
        self.k_th = alpha2_lon/b_theta
        self.k_thdot = alpha1_lon/b_theta
        self.ki_lon = -alpha3_lon/b_theta
        alpha1_lat = pi_psi + 2*wn_phi*zeta_phi + 2*wn_psi*zeta_psi
        alpha2_lat = 2*pi_psi*wn_phi*zeta_phi + 2*pi_psi*wn_psi*zeta_psi + wn_phi**2 + 4*wn_phi*wn_psi*zeta_phi*zeta_psi + wn_psi**2
        alpha3_lat = pi_psi*wn_phi**2 + 4*pi_psi*wn_phi*wn_psi*zeta_phi*zeta_psi + pi_psi*wn_psi**2 + 2*wn_phi**2*wn_psi*zeta_psi + 2*wn_phi*wn_psi**2*zeta_phi
        alpha4_lat = 2*pi_psi*wn_phi**2*wn_psi*zeta_psi + 2*pi_psi*wn_phi*wn_psi**2*zeta_phi + wn_phi**2*wn_psi**2
        alpha5_lat = pi_psi*wn_phi**2*wn_psi**2
        b1 = 1/P.J1x
        a1 = P.ellT*Fe/(JT+P.J1z)
        self.k_phi = alpha2_lat/b1
        self.k_psi = alpha4_lat/(a1*b1) #?why a1*b1? Ans: my version of the a1 term included 1/J1x, his doesn't so he multiplies
        self.k_phidot = alpha1_lat/b1
        self.k_psidot = alpha3_lat/(a1*b1)
        self.ki_lat = -alpha5_lat/(a1*b1)
        # print gains to terminal
        print('K_lon: [', self.k_th, ',', self.k_thdot, ']')
        print('ki_lon: ', self.ki_lon)         
        print('K_lat: [', self.k_phi, ',', self.k_psi, ',', self.k_phidot, ',', self.k_psidot, ']')
        print('ki_lat: ', self.ki_lat)        
        #--------------------------------------------------
        # saturation limits
        theta_max = 30.0 * np.pi / 180.0  # Max theta, rads
        #--------------------------------------------------
        self.Ts = P.Ts
        self.sigma = 0.05  # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - self.Ts) / (2 * self.sigma + self.Ts)
        self.phi_d1 = 0.
        self.phi_dot = 0.
        self.theta_d1 = 0.
        self.theta_dot = 0.
        self.psi_d1 = 0.
        self.psi_dot = 0.        
        # variables to implement integrator
        self.integrator_th = 0.0  
        self.error_th_d1 = 0.0  
        self.integrator_psi = 0.0  
        self.error_psi_d1 = 0.0 

    def update(self, r, y):
        theta_ref = r[0][0]
        psi_ref = r[1][0]
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]
        force_equilibrium = (P.m1*P.ell1 + P.m2*P.ell2)*P.g*np.cos(theta) / P.ellT
        # update differentiators
        self.phi_dot = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.phi_dot   + 2.0 / (2.0*self.sigma + self.Ts) * (phi - self.phi_d1)
        self.phi_d1 = phi
        self.theta_dot = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.theta_dot + 2.0 / (2.0*self.sigma + self.Ts) * (theta - self.theta_d1)
        self.theta_d1 = theta
        self.psi_dot = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.psi_dot   + 2.0 / (2.0*self.sigma + self.Ts) * (psi - self.psi_d1)  
        self.psi_d1 = psi
        # integrate error
        error_th = theta_ref - theta
        error_psi = psi_ref - psi
        self.integrator_th = self.integrator_th + (self.Ts / 2.0) * (error_th + self.error_th_d1)
        self.integrator_psi = self.integrator_psi + (self.Ts / 2.0) * (error_psi + self.error_psi_d1)
        self.error_th_d1 = error_th
        self.error_psi_d1 = error_psi

        # longitudinal control
        force_unsat = force_equilibrium - self.k_th *theta - self.k_thdot *self.theta_dot - self.ki_lon *self.integrator_th
        force = saturate(force_unsat, -P.force_max, P.force_max)
        # lateral control
        torque_unsat = -self.k_phi *phi - self.k_psi *psi - self.k_phidot *self.phi_dot - self.k_psidot *self.psi_dot - self.ki_lat *self.integrator_psi
        torque = saturate(torque_unsat, -P.torque_max, P.torque_max)
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)
        return pwm, np.array([[0], [theta_ref], [psi_ref]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u
