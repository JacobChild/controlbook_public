import numpy as np
import control as cnt
import hummingbirdParam as P
import sympy as sym

#Work for prb 12
#1
#longitudinal symbolic equations
eqn1Lon = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_th')*sym.Symbol('wn_th')*sym.Symbol('s') + sym.Symbol('wn_th')**2
eqn2Lon = sym.Symbol('s') + sym.Symbol('pi_th')
eqn3Lon = eqn1Lon* eqn2Lon
expandedLon = sym.expand(eqn3Lon)
coefLon = sym.collect(expandedLon, sym.Symbol('s'))
print(coefLon)
#lateral symbolic equations
eqn1Lat = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_phi')*sym.Symbol('wn_phi')*sym.Symbol('s') + sym.Symbol('wn_phi')**2
eqn2Lat = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_psi')*sym.Symbol('wn_psi')*sym.Symbol('s') + sym.Symbol('wn_psi')**2
eqn3Lat = sym.Symbol('s') + sym.Symbol('pi_psi')
eqn4Lat = eqn1Lat* eqn2Lat* eqn3Lat
expandedLat = sym.expand(eqn4Lat)
coefLat = sym.collect(expandedLat, sym.Symbol('s'))
print(coefLat)

#2
JT = P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J2z + P.m3 * (P.ell3x**2 + P.ell3y**2)
Fe = (P.m1*P.ell1 + P.m2*P.ell2)*P.g / P.ellT
b_theta = P.ellT/(P.m1 * P.ell1**2 + P.m2 * P.ell2**2 + P.J1y + P.J2y)
ALon = np.array([[0.0,1.0],
                 [0.0,0.0]])
BLon = np.array([[0.0],[b_theta]])
CrLon = np.array([[1.0,0.0]])
ALatTerm1 = P.ellT * Fe / (JT + P.J1z)
ALat = np.array([[0.0,0.0,1.0,0.0],
                 [0.0,0.0,0.0,1.0],
                 [0.0,0.0,0.0,0.0],
                 [ALatTerm1,0.0,0.0,0.0]])
BLat = np.array([[0.0],
                 [0.0],
                 [1.0/(P.J1x)],
                 [0.0]])
CrLat = np.array([[0.0,1.0,0.0,0.0]])

ALonAug = np.vstack((np.hstack((ALon, np.zeros((np.size(ALon,1),1)))),
                               np.hstack((-CrLon, np.array([[0.0]]))) ))
BLonAug = np.vstack((BLon, 0.0))

ALatAug = np.vstack((np.hstack((ALat, np.zeros((np.size(ALat,1),1)))),
                     np.hstack((-CrLat, np.array([[0.0]]))) ))
BLatAug = np.vstack((BLat, 0.0))

#3 Create the controlability matrix
CLoncnt = cnt.ctrb(ALonAug, BLonAug)
print("\n")
print(CLoncnt)
print("\n")
CLatcnt = cnt.ctrb(ALatAug, BLatAug)
print(CLatcnt)
#inverse them
CLoncntInv = np.linalg.inv(CLoncnt)
CLatcntInv = np.linalg.inv(CLatcnt)
print("\n")
print(CLoncntInv)
print("\n")
print(CLatcntInv)

#4 check the determinants of the controlability matrices and see if they are s^3 and s^5 respectively
#? is that the same as checking rank?
print("\n")
if np.linalg.matrix_rank(CLoncnt) == 3 and np.linalg.matrix_rank(CLatcnt) == 5:
    print("The system is controllable")
else:
    print("The system is not controllable")