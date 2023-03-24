import numpy as np
import VTOLParam as P

# set up the VTOL dynamics class
class VTOLDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([[P.z0], [P.h0], [P.theta0], [P.zdot0], [P.hdot0], [P.thetadot0]])  # initial position and velocity for z, h, and theta
        #Masses of the VTOL, kg
        self.mc = P.mc * (1. + alpha * (2. * np.random.rand() - 1.)) # mass of the cabin
        self.mr = P.mr * (1. + alpha * (2. * np.random.rand() - 1.)) # mass of the rotor
        self.ml = P.ml * (1. + alpha * (2. * np.random.rand() - 1.)) # mass of the payload
        # Other VTOL Parameters
        self.Jc = P.Jc * (1. + alpha * (2. * np.random.rand() - 1.)) # moment of inertia of the cabin
        self.d = P.d * (1. + alpha * (2. * np.random.rand() - 1.)) # distance from the rotor to the cabin
        # Environment Parameters
        self.mu = P.mu #* (1. + alpha * (2. * np.random.rand() - 1.)) # dynamic viscosity of air
        self.g = P.g # gravity
        self.F_wind = P.F_wind * (1. + alpha * (2. * np.random.rand() - 1.)) # wind disturbance force
        # Sample rate at which the dynamics are propagated
        self.Ts = P.Ts
        self.fmax = P.fmax
        
    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        #u = self.saturate(u, self.fmax) #calls the saturate function
        self.rk4_step(u) #propagate the state by one time step
        y = self.h()  # return the corresponding output
        return y
    
    def f(self, state, u):
        # Return xdot = f(x,u), the system state update equations
        #convert the VTOL system to state space
        fr = u[0][0]
        fl = u[1][0]
        z = state[0][0]
        h = state[1][0]
        theta = state[2][0]
        zdot = state[3][0]
        hdot = state[4][0]
        thetadot = state[5][0]
        zddot = ( -1.0*(fr + fl)*np.sin(theta) - self.mu*zdot) / (self.mc + 2.0*self.mr) #? for the future - self.mu*zdot
        hddot = ( -1.0*(self.mc + 2.0*self.mr)*self.g + (fr + fl)*np.cos(theta)) / (self.mc + 2.0*self.mr)
        thetaddot = ( self.d * (fr - fl) ) / (self.Jc + 2.0*self.mr*self.d**2)
        xdot = np.array([[zdot], [hdot], [thetadot], [zddot], [hddot], [thetaddot]])
        return xdot
    
    def h(self):
        #return the output equations
        #could als use input u if needed
        z = self.state[0][0]
        h = self.state[1][0]
        theta = self.state[2][0]
        y = np.array([[z], [h], [theta]])
        return y
    
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4) 
        
    
    def saturate(self, u, limit):
        # Saturate the input u
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u