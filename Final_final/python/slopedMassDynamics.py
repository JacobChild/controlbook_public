import numpy as np
import slopedMassParam as P

class slopedMassDynamics:
    def __init__(self, alpha=0.0):
        self.state = np.array([
            [0.0],
            [0.0]
        ])
        self.g = 9.8
        self.theta = 45 * np.pi / 180 * (1.+alpha*(2.*np.random.rand()-1.))
        self.m = 0.5 * (1.+alpha*(2.*np.random.rand()-1.))
        self.k1 = 0.05 * (1.+alpha*(2.*np.random.rand()-1.))
        self.k2 = 0.02 * (1.+alpha*(2.*np.random.rand()-1.))
        self.b = 0.1 * (1.+alpha*(2.*np.random.rand()-1.))
        self.Ts = 0.01
        self.force_limit = P.F_max

    def update(self, u):
        u = saturate(u, self.force_limit)
        self.rk4_step(u)
        y = self.h()
        return y

    def f(self, state, F):
        z = state[0][0]
        zdot = state[1][0]
        xdot = np.array([
            [zdot],
            [-(self.k1 / self.m) * z \
             - (self.k2 / self.m) * z**3 \
             - (self.b / self.m) * zdot \
             + self.g * np.sin(self.theta) \
             + saturate(F, self.force_limit) / self.m]
        ])
        return xdot

    def h(self):
        z = self.state[0][0]
        return z

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u