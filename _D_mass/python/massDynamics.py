import numpy as np
import massParam as P

class massDynamics:
    def __init__(self, alpha = 0.0):
        # Initial state conditions
        self.state = np.array([[P.z0], [P.zdot0]]) #initial position and velocity
        #Mass of the block, kg
        self.m = P.m * (1+2*alpha*np.random.rand()-alpha)
        #Spring constant, N/m
        self.k = P.k * (1+2*alpha*np.random.rand()-alpha)
        #Damping coefficient, Ns/m
        self.b = P.b * (1+2*alpha*np.random.rand()-alpha)
        #Sample rate at which the dynamics are propagated
        self.Ts = P.Ts
        self.force_limit = P.F_max
        
    def update(self, u):
        #This is the external method that takes the input u at time t and returns the output y at time t.
        #Saturate the force
        u = self.saturate(u, P.F_max) #calls the saturate function and updates u if needed #! make the saturate function
        self.rk4_step(u) #propagate the state by one time sample
        y = self.h() #return the corresponding output
        return y
    
    def f(self, state, u): #This is my Equation of motion
        # Return xdot = f(x,u), the system state update equations
        # convert the mass spring damping system to state space
        z = state[0][0]
        zdot = state[1][0]
        zddot = (1.0/self.m)*(u - self.b*zdot - self.k*z)
        xdot = np.array([[zdot], [zddot]])
        return xdot
    
    def h(self):
        #return the output equations
        #could also use input u if needed
        z = self.state[0][0]
        y = np.array([[z]])
        return y
    
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)    
        
    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u) #if the input torque is larger than the limit, then the input torque is set to the limit
        return u