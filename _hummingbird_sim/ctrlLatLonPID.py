import numpy as np
import hummingbirdParam as P


class ctrlLatLonPID:
    def __init__(self):
        # tuning parameters
        tr_pitch = 0.650 # rise time for pitch
        zeta_pitch = .707 # damping ratio for pitch
        self.ki_pitch = 1.1 # integrator gain
        self.ki_yaw = 0.0001
        tr_phi = .3
        zeta_phi = .7
        tr_psi = tr_phi * 10.0 # rise time for yaw (outer loop) needs to be faster than inner
        zeta_psi = .7
        #other needed terms
        JT = P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J2z + P.m3 * (P.ell3x**2 + P.ell3y**2)
        Fe = (P.m1*P.ell1 + P.m2*P.ell2)*P.g / P.ellT
        # gain calculation
        b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
        b_phi = 1.0/P.J1x
        b_psi = P.ellT * Fe / (JT + P.J1z)
        #print('b_theta: ', b_theta)
        wn_pitch = 2.2/tr_pitch # natural frequency for pitch
        wn_phi = 2.2/tr_phi
        wn_psi = 2.2/tr_psi
        self.kp_pitch = wn_pitch**2/b_theta
        self.kd_pitch = 2.0*zeta_pitch*wn_pitch/b_theta
        self.kp_phi = wn_phi**2/b_phi
        self.kd_phi = 2.0*zeta_phi*wn_phi/b_phi
        self.kp_psi = wn_psi**2/b_psi
        self.kd_psi = 2.0*zeta_psi*wn_psi/b_psi
        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('ki_pitch: ', self.ki_pitch)
        print('kd_pitch: ', self.kd_pitch) 
        print('kp_phi: ', self.kp_phi)
        print('kd_phi: ', self.kd_phi)
        print('kp_psi: ', self.kp_psi)
        print('kd_psi: ', self.kd_psi)
        
        # sample rate of the controller
        self.Ts = P.Ts
        # dirty derivative parameters
        self.sigma = 0.005  # cutoff freq for dirty derivative
        self.beta = (2 * self.sigma - self.Ts) / (2 * self.sigma + self.Ts)
        # delayed variables
        self.phi_d1 = 0.
        self.theta_d1 = 0.
        self.psi_d1 = 0.
        self.phi_dot = 0.
        self.theta_dot = 0.
        self.psi_dot = 0.
        self.integrator_theta = 0.
        self.integrator_psi = 0.
        self.error_theta_d1 = 0.  # pitch error delayed by 1
        self.error_psi_d1 = 0.  # roll error delayed by 1

    def update(self, r, y):
        #print('r: ', r)
        #print('y: ', y)
        theta_ref = r[0][0]
        psi_ref = r[1][0]
    
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]
        
        force_fl = (P.m1*P.ell1 + P.m2*P.ell2)*P.g*np.cos(theta) / P.ellT
        # compute errors
        error_theta = theta_ref - theta
        error_psi = psi_ref - psi
        # update differentiators
        self.theta_dot = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.theta_dot + 2.0 / (2.0*self.sigma + self.Ts) * (theta - self.theta_d1)
        self.psi_dot   = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.psi_dot   + 2.0 / (2.0*self.sigma + self.Ts) * (psi - self.psi_d1)
        self.phi_dot   = (2.0*self.sigma - self.Ts) / (2.0*self.sigma + self.Ts) * self.phi_dot   + 2.0 / (2.0*self.sigma + self.Ts) * (phi - self.phi_d1)
        # update integrators
        #integrator anti windup
        #print("thetadot: ", self.theta_dot)
        if np.abs(self.theta_dot) < 0.01:
            #print("integrator on")
            self.integrator_theta = self.integrator_theta + (self.Ts / 2.0) * (error_theta + self.error_theta_d1)
            
        # pitch control
        force_unsat = self.kp_pitch * error_theta + self.ki_pitch * self.integrator_theta - self.kd_pitch * self.theta_dot + force_fl
        #force_unsat = force_fl
        force = saturate(force_unsat, -P.force_max, P.force_max)
        
        # yaw control
        #integrator and anti windup
        if np.abs(self.psi_dot) < 0.01:
            #print("yaw integrator on")
            self.integrator_psi = self.integrator_psi + (self.Ts / 2.0) * (error_psi + self.error_psi_d1)
            
        phi_ref = self.kp_psi * error_psi + self.ki_yaw * self.integrator_psi - self.kd_psi * self.psi_dot
        error_phi = phi_ref - phi
            
        T_phi = self.kp_phi * error_phi - self.kd_phi * self.phi_dot
        torque = T_phi 
        
        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # u_right          
        pwm = saturate(pwm, 0, 0.7)
        # update all delayed variables
        self.theta_d1 = theta
        self.phi_d1 = phi
        self.psi_d1 = psi
        self.error_theta_d1 = error_theta

        # return pwm plus reference signals
        return pwm, np.array([[phi_ref], [theta_ref], [psi_ref]]), np.array([[force], [torque]])


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




