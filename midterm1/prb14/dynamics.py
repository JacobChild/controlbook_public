import numpy as np 
import param as P


class Dynamics:
    def __init__(self):
        self.state = np.array([
            [P.y0],      
            [P.ydot0]    
        ])  
        self.Ts = P.Ts  

    def update(self, u):
        self.rk4_step(u)  
        y = self.h()  
        return y

    def f(self, state, u):
        y = state[0][0]
        ydot = state[1][0]
        yddot = 3 * u + 2 * ydot
        xdot = np.array([[ydot], [yddot]])
        return xdot

    def h(self):
        # return the output equations
        # could also use input u if needed
        y = self.state[0][0]
        y = np.array([[y]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
