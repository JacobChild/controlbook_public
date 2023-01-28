import numpy as np
import blockbeamParam as P


# set up the block beam dynamics class
class blockBeamDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([[P.z0], [P.theta0], [P.zdot0], [P.thetadot0]])  # initial position and velocity for z and theta
        # Mass of the block, kg
        self.mbl = P.m1 * (1. + alpha * (2. * np.random.rand() - 1.)) # mass of the block
        self.mbe = P.m2 * (1. + alpha * (2. * np.random.rand() - 1.)) # mass of the beam
        # Length of the beam, m
        self.ell = P.length * (1. + alpha * (2. * np.random.rand() - 1.)) # length of the beam
        # the gravity constant is well known, so we don't change it.
        self.g = P.g
        # sample rate at which the dynamics are propagated
        self.Ts = P.Ts
        self.force_limit = P.Fmax
        
    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = self.saturate(u, self.force_limit) #calls the saturate function
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, u): #This is my Equation of motion
        # Return xdot = f(x,u), the system state update equations
        # convert the block beam system to state space
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]
        zddot = z*thetadot**2 - self.g*np.sin(theta)
        thetaddot = (-2*self.mbl*z*zdot*thetadot - self.mbl*self.g*z*np.cos(theta) - self.mbe*self.g*self.ell/2.0 * np.cos(theta) + u*self.ell*np.cos(theta)) / (self.mbe*self.ell**2/3 + self.mbl*z**2)
        xdot = np.array([[zdot], [thetadot], [zddot], [thetaddot]])
        return xdot
    
    def h(self):
        #return the output equations
        #could also use input u if needed
        z = self.state[0]
        theta = self.state[2]
        y = np.array([z, theta])
        return y
    
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4) 
        
    def saturate(self, u, limit):
        # saturate the input u
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u