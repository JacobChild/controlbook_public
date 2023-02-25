import numpy as np

class Dynamics():
    def __init__(self):
        # Initial state conditions
        self.J = 0.1
        self.r = 2.0
        self.bigR = 3.2
        self.m = 3.3
        self.b = .25
        self.Force_max = 100.0
        self.phi0 = 0.0
        self.z0 = 0.0
        self.phidot0 = 0.0
        self.zdot0 = 0.0
        self.state = np.array([self.phi0], [self.z0], [self.phidot0], [self.zdot0])
        
    def update(self, u):
         # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = self.saturate(u, self.Force_max)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
        
    def f(self, state, u):
        phi = state[0][0]
        z = state[1][0]
        phidot = state[2][0]
        zdot = state[3][0]
        phiddot = (1/(1-self.r**2/(self.J * self.m)))*(self.bigR*u/self.J + z*np.sin(phi) - self.r/(self.m*self.J)*(-self.b*zdot - np.cos(phi))/(self.m*self.J))
        zddot = (1/(1-self.r**2/self.J))*(-self.b*zdot - np.cos(phi)- self.r*(self.bigR*u+z*np.sin(phi))/self.J)
        xdot = np.array([[phidot], [zdot], [phiddot], [zddot]])
        return xdot
        
    def h(self):
        # return the output equations
        # could also use input u if needed
        phi = self.state[0][0]
        z = self.state[1][0]
        y = np.array([[phi], [z]])
        return y
        
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state = self.state + self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
            
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u