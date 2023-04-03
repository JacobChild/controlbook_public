import numpy as np
#import sympy as sym
import sympy as sym
import control as cnt

#Problem 2
# create a 4x4 symbolic alphabet matrix
A = sym.Matrix([[sym.Symbol('s-1'), sym.Symbol('-b'), sym.Symbol('c'), sym.Symbol('d')],
                [-1, sym.Symbol('s'), 0, 0],
                [0, -1, sym.Symbol('s'), 0],
                [0, 0, -1, sym.Symbol('s')]])
#print(A)
#A2 = sym.Matrix.inv(A)
#print(A2)

#problem 5
# check for linear independence
A = np.array([[1,2,3,4],
              [5,6,7,8],
              [3*2,6*2,9*2,12*2],
              [4*3,8*3,12*3,16*3]])
#check for controllability
if np.linalg.matrix_rank(A) == np.size(A,1):
    print('the system is controllable')
else: 
    print('the system is not controllable')
    
#Problem 8
A = np.array([[-1.0, 4.0],
              [1.0, 0.0]])
B = np.array([[1.0], [0.0]])
C = np.array([1.0, -1.0])

K = cnt.place(A, B, [-1.0 + 2.j, -1.0 - 2.j])
print('K: ', K)

#Problem 9
K = np.array([[2.0, 5.0]])
kr = -1.0 / (C @ np.linalg.inv(A - B @ K) @ B)
print('kr: ', kr)

#Problem 10 #? I am unsure on this one
A1 = np.array([[-5.0, 0.0],
               [-1.0, 0.0]])
B1 = np.array([[1.0], [0.0]])
CC = cnt.ctrb(A1, B1)
print(CC)
if np.linalg.matrix_rank(CC) == np.size(A1,1):
    print('the system is controllable, prb 10')
else: 
    print('the system is not controllable, problem 10')
alpha = sym.Matrix([[sym.Symbol('p')+1, sym.Symbol('p')]])
littleaBigA = sym.Matrix([[5, 0]])
curlyA = sym.Matrix([[1, 5], [0, 1]])
CCsym = sym.Matrix([[1, -5], [0, -1]])
K1 = (alpha - littleaBigA) @ sym.Matrix.inv(curlyA) @ sym.Matrix.inv(CCsym)
print("K1: ", K1)
