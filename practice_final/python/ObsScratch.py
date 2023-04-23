import numpy as np
import sympy as sym

Amat = sym.Matrix([[sym.Symbol('aA'), sym.Symbol('aB')],
               [sym.Symbol('aC'), sym.Symbol('aD')]])
Bmat = sym.Matrix([[sym.Symbol('bA')],
               [sym.Symbol('bB')]])
Cmat = sym.Matrix([[sym.Symbol('cA'), sym.Symbol('cB')]])
A1 = np.vstack((np.hstack((Amat, np.zeros((np.size(Amat, 1),1)))),
                        np.hstack((-Cmat, np.array([[0.0]]))) ))
B1 = np.vstack((Bmat, 0.0))
print("A1: ",A1)
print("\n")
print("B1",B1)