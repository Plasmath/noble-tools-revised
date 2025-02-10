#vertices.py
#Code that handles the generation of orbits and orbit types.

import numpy as np
import sympy as sp
import config

phi = ( 1 + np.sqrt(5) ) / 2
root2 = np.sqrt(2)

#Generating reflection matrices for the *332 reflection group
r1tet = np.array([[ 0, 1, 0],
                  [ 1, 0, 0],
                  [ 0, 0, 1]])
r2tet = np.array([[ 0,-1, 0],
                  [-1, 0, 0],
                  [ 0, 0, 1]])
r3tet = np.array([[ 1, 0, 0],
                  [ 0, 0, 1],
                  [ 0, 1, 0]])

#Generating reflection matrices for the *432 reflection group
r1oct = np.array([[-1, 0, 0],
                  [ 0, 1, 0],
                  [ 0, 0, 1]])
r2oct = r1tet
r3oct = r3tet

#Generating reflection matrices for the *532 reflection group
r1ico = r1oct
r2ico = np.array([[ 1, 0, 0],
                  [ 0,-1, 0],
                  [ 0, 0, 1]])
r3ico = 1/2 * np.array([[1-phi,  -phi,     1],
                        [ -phi,     1, phi-1],
                        [    1, phi-1,   phi]])

#basis vectors for *332 reflection group
v1tet = np.array([-root2/2, root2/2, root2/2])
v2tet = np.array([ root2/2, root2/2, root2/2])
v3tet = np.array([ 0, 0, root2])

#basis vectors for *432 reflection group
v1oct = np.array([1,1,1])
v2oct = np.array([0,root2,root2])
v3oct = np.array([0,0,root2])

#basis vectors for *532 reflection group
v1ico = np.array([1,0,phi+1])
v2ico = np.array([0,1,phi])
v3ico = np.array([0,0,2*phi])

def InitializeSymbolic(): #initialize basis vectors and matrices for symbolic mode
    global phi
    global root2
    global r3ico
    global v1tet
    global v2tet
    global v3tet
    global v2oct
    global v3oct
    global v1ico
    global v2ico
    global v3ico

    phi = ( 1 + sp.sqrt(5) ) / 2
    root2 = sp.sqrt(2)
    
    r3ico = sp.Rational(1,2) * np.array([[1-phi,  -phi,     1],
                                         [ -phi,     1, phi-1],
                                         [    1, phi-1,   phi]])
    v1tet = np.array([-root2/2, root2/2, root2/2])
    v2tet = np.array([ root2/2, root2/2, root2/2])
    v3tet = np.array([ 0, 0, root2])
    v2oct = np.array([0,root2,root2])
    v3oct = np.array([0,0,root2])
    v1ico = np.array([1,0,phi+1])
    v2ico = np.array([0,1,phi])
    v3ico = np.array([0,0,2*phi])

def RoundArray(array):
    return array.round(decimals = config.numericPrecision)

def IsNew(p,pts):
    if config.useSymbolic:
        return not any(np.all(p == q) for q in pts)
    else:
        return not np.any(np.all(RoundArray(p) == [RoundArray(k) for k in pts], axis=1))

def GenerateFromGroup(initial,r1,r2,r3,N):
    if config.useSymbolic:
        pts = [[sp.expand(k) for k in initial]]
    else:
        pts = [initial]
    
    for i in range(N):
        for p in pts:
            if config.useSymbolic:
                p1 = [sp.expand(k) for k in r1.dot(p)]
                p2 = [sp.expand(k) for k in r2.dot(p)]
                p3 = [sp.expand(k) for k in r3.dot(p)]
            else:
                p1 = r1.dot(p)
                p2 = r2.dot(p)
                p3 = r3.dot(p)
            
            if IsNew(p1,pts): pts.append(p1)
            if IsNew(p2,pts): pts.append(p2)
            if IsNew(p3,pts): pts.append(p3)
    return pts

def GenerateVStar332(a,b,c): #Generate orbit of V_*332 with parameters a,b,c.
    initial = a*v1tet + b*v2tet + c*v3tet
    return GenerateFromGroup(initial, r1tet, r2tet, r3tet, 1)

def GenerateV332(a,b,c): #Generate orbit of V_332 with parameters a,b,c
    initial = a*v1tet + b*v2tet + c*v3tet
    return GenerateFromGroup(initial, r1tet.dot(r2tet), r1tet.dot(r3tet), r2tet.dot(r3tet), 1)

def GenerateV3Star2(a,b,c): #Generate orbit of V_3*2 with parameters a,b,c
    initial = a*v1oct + b*v2oct + c*v3oct
    return GenerateFromGroup(initial, r1oct, r2oct.dot(r1oct.dot(r2oct)), r3oct.dot(r2oct), 1)

def GenerateVStar432(a,b,c): #Generate orbit of V_*432 with parameters a,b,c.
    initial = a*v1oct + b*v2oct + c*v3oct
    return GenerateFromGroup(initial, r1oct, r2oct, r3oct, 1)

def GenerateV432(a,b,c): #Generate orbit of V_432 with parameters a,b,c
    initial = a*v1oct + b*v2oct + c*v3oct
    return GenerateFromGroup(initial, r1oct.dot(r2oct), r1oct.dot(r3oct), r2oct.dot(r3oct), 1)

def GenerateVStar532(a,b,c): #Generate orbit of V_*532 with parameters a,b,c.
    initial = a*v1ico + b*v2ico + c*v3ico
    return GenerateFromGroup(initial, r1ico, r2ico, r3ico, 1)

def GenerateV532(a,b,c): #Generate orbit of V_532 with parameters a,b,c
    initial = a*v1ico + b*v2ico + c*v3ico
    return GenerateFromGroup(initial, r1ico.dot(r2ico), r1ico.dot(r3ico), r2ico.dot(r3ico), 1)