import numpy as np
#import param as P

class ctrlPD:
    def __init__(self):
        zeta = .707 #damping ratio
        tr = 0.50 #rise time
        wn = 2.2/tr #omegan, natural frequency
        b0 = 3.0 #from transfer function
        a0 = 0.0 #from transfer function
        a1 = -2.0 #from transfer function
        self.kp = (wn**2.0 - a0) / b0
        self.kd = (2.0*zeta*wn - a1) / b0
        print('kp = ', self.kp)
        print('kd = ', self.kd)

    def update(self, ref, state):
        y = state[0][0]
        yerr = ref - y
        ydot = state[1][0]
        u = self.kp*yerr - self.kd*ydot
        return u
