import sympy as sym
from sympy import Inverse
import numpy as np
import control as cnt

#Work for prb 12
#1
#longitudinal symbolic equations
eqn1Lon = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_th')*sym.Symbol('wn_th')*sym.Symbol('s') + sym.Symbol('wn_th')**2
eqn2Lon = sym.Symbol('s') + sym.Symbol('pi_th')
eqn3Lon = eqn1Lon* eqn2Lon
expandedLon = sym.expand(eqn3Lon)
coefLon = sym.collect(expandedLon, sym.Symbol('s'))
print("\n")
print(coefLon)
print("\n")
#lateral symbolic equations
eqn1Lat = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_phi')*sym.Symbol('wn_phi')*sym.Symbol('s') + sym.Symbol('wn_phi')**2
eqn2Lat = sym.Symbol('s')**2 + 2*sym.Symbol('zeta_psi')*sym.Symbol('wn_psi')*sym.Symbol('s') + sym.Symbol('wn_psi')**2
eqn3Lat = sym.Symbol('s') + sym.Symbol('pi_psi')
eqn4Lat = eqn1Lat* eqn2Lat* eqn3Lat
expandedLat = sym.expand(eqn4Lat)
coefLat = sym.collect(expandedLat, sym.Symbol('s'))
print(coefLat)

JT = sym.Symbol('JT')
Fe = sym.Symbol('Fe')
b_theta = sym.Symbol('b_theta')
ellT = sym.Symbol('ellT')
J1z = sym.Symbol('J1z')
J1x = sym.Symbol('J1x')
s = sym.Symbol('s')
ALon = sym.Matrix([[0.0,1.0],
                 [0.0,0.0]])
BLon = sym.Matrix([[0.0],[b_theta]])
CrLon = sym.Matrix([[1.0,0.0]])
ALatTerm1 = Fe * ellT / (JT + J1z)
ALat = sym.Matrix([[0.0,0.0,1.0,0.0],
                 [0.0,0.0,0.0,1.0],
                 [0.0,0.0,0.0,0.0],
                 [ALatTerm1,0.0,0.0,0.0]])
BLat = sym.Matrix([[0.0],
                   [0.0],
                   [1.0/J1x],
                   [0.0]])
CrLat = sym.Matrix([[0.0,1.0,0.0,0.0]])

ALonAug = np.vstack((np.hstack((ALon, np.zeros((np.size(ALon,1),1)))),
                               np.hstack((-CrLon, np.array([[0.0]]))) ))
BLonAug = np.vstack((BLon, 0.0))

ALatAug = np.vstack((np.hstack((ALat, np.zeros((np.size(ALat,1),1)))),
                     np.hstack((-CrLat, np.array([[0.0]]))) ))
BLatAug = np.vstack((BLat, 0.0))

#print the  shape of each of the above terms
print(sym.shape(ALonAug), "  ", sym.shape(BLonAug), "  ", sym.shape(CrLon))
print("\n")
print(ALonAug@BLonAug)


CLoncnt = sym.Matrix([BLonAug, ALonAug @ BLonAug, (ALonAug@ALonAug) @ BLonAug])
CLatcnt = sym.Matrix([BLatAug, ALatAug @ BLatAug, (ALatAug@ALatAug) @ BLatAug, (ALatAug@ALatAug@ALatAug) @BLatAug, (ALatAug@ ALatAug@ALatAug@ALatAug) @ BLatAug])
print(CLoncnt)
print("\n")
print(CLatcnt)
print("\n")

LonIMatrix = sym.Matrix([[1,0,0],
                         [0,1,0],
                         [0,0,1]])
LatIMatrix =  sym.Matrix ([[1,0,0,0,0],
                           [0,1,0,0,0],
                           [0,0,1,0,0],
                           [0,0,0,1,0],
                           [0,0,0,0,1]])
LonRank = sym.Matrix.det(s*LonIMatrix - ALonAug)
LatRank = sym.Matrix.det(s*LatIMatrix - ALatAug)
print(LonRank)
print("\n")
print(LatRank)

aLon = sym.Matrix([[0],[0],[0]])
aLat = sym.Matrix([[0],[0],[0],[0],[0]])
alphaLon = sym.Matrix([[sym.Symbol('alpha1Lon')], [sym.Symbol('alpha2Lon')], [sym.Symbol('alpha3Lon')]])
alphaLat = sym.Matrix([[sym.Symbol('alpha1Lat')], [sym.Symbol('alpha2Lat')], [sym.Symbol('alpha3Lat')], [sym.Symbol('alpha4Lat')], [sym.Symbol('alpha5Lat')]])
CurlyAinvLon = LonIMatrix
CurlyAinvLat = LatIMatrix
print("\n")

#hard code the inverted matrix
A = sym.Symbol('A')
B = sym.Symbol('B')
CTestLon = sym.Matrix([[0,b_theta,0],
                       [b_theta,0,0],
                       [0,0,-b_theta]])
CTestLat = sym.Matrix([[0,0,A,0,0],
                       [A,0,0,0,0],
                       [0,0,0,B,0],
                       [0,B,0,0,0],
                       [0,0,0,0,-B]])



KLon = (alphaLon - aLon).T @ CurlyAinvLon @ sym.Matrix.inv(CTestLon)
KLat = (alphaLat - aLat).T @ CurlyAinvLat @ sym.Matrix.inv(CTestLat)
print("\n")
print(KLon)
print("\n")
print(KLat)