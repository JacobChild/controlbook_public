import massParam as P
import numpy as np
from digitalFilter import digitalFilter
import loopShaping as L

class ctrlLoopshape:
    def __init__(self, method="state_space"):
        self.Ts = P.Ts
        if method == "state_space":
            self.x_C = np.zeros((L.C_ss.A.shape[0], 1))
            self.x_F = np.zeros((L.F_ss.A.shape[0], 1))
            self.A_C = L.C_ss.A
            self.B_C = L.C_ss.B
            self.C_C = L.C_ss.C
            self.D_C = L.C_ss.D
            self.A_F = L.F_ss.A
            self.B_F = L.F_ss.B
            self.C_F = L.F_ss.C
            self.D_F = L.F_ss.D
            self.N = 10  #number of Euler integration steps for each sample
        elif method == "digital_filter":
            self.prefilter = digitalFilter(L.F.num, L.F.den, P.Ts)
            self.control = digitalFilter(L.C.num, L.C.den, P.Ts)
        self.method = method

    def update(self, z_r, y):
        zCurrent = y[0][0]
        # prefilter
        if self.method == "state_space":
            self.updatePrefilterState(z_r)
            z_r_filtered = (self.C_F @ self.x_F + self.D_F * z_r)[0,0]
        elif self.method == "digital_filter":
            z_r_filtered = self.prefilter.update(z_r)
        # filtered error signal
        error = z_r_filtered - zCurrent
        # update controller
        if self.method == "state_space":
            self.updateControlState(error)
            F_tilde = (self.C_C @ self.x_C + self.D_C * error)[0,0]
        elif self.method == "digital_filter":
            F_tilde = self.control.update(error)
        # compute feedback linearization torque F_tilde
        #linearized force?
        # compute total torque
        Force = saturate(F_tilde, P.F_max)
        return Force.item(0)

    def updatePrefilterState(self, z_r):
        for i in range(0, self.N):
            self.x_F = self.x_F + (self.Ts/self.N)*(
                self.A_F @ self.x_F + self.B_F * z_r
            )

    def updateControlState(self, error_out):
        for i in range(0, self.N):
            self.x_C = self.x_C + (self.Ts/self.N)*(
                self.A_C @ self.x_C + self.B_C * error_out
            )


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
