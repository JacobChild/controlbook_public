import control as cnt
import numpy as np

A = np.array([[0.0, 1.0],
              [2.0, 3.0]])
B = np.array([[0.0], [4.0]])
CC = cnt.ctrb(A, B)
CCTest = np.array([B, A@B]) #! this has the same numbers, but a bit of a different format than above
Cr = np.array([5.0, 0.0])

#check if they are the same
#print(CC.T)
#print(CCTest)

#check for controllability
if np.linalg.matrix_rank(CC) == np.size(A,1):
    print('the system is controllable')
else:
    print('the system is not controllable')
    
asubA = np.array([-3.0, -2.0])
curlyAsubA = np.array([[1.0, -3.0], [0.0, 1.0]])
alpha = np.array([2.0, 2.0])
K = (alpha - asubA) @ np.linalg.inv(curlyAsubA) @ np.linalg.inv(CC)
print('B: ', B)
Kr = -1.0 / (Cr @ np.linalg.inv(A - B @ K) @ B)
print('K: ', K)
print('Kr: ', Kr)

#code way -
KCode = cnt.place(A, B, [-1.0, 1.0])
print('KCode: ', KCode)