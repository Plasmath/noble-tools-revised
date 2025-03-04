#Code that constructs the symbolic forms of the orbit types.
#This is only needed for the orbit types with one or more degrees of freedom
#in order to obtain exact minimal polynomials for their parameters.
import vertices
import config
import sympy as sp
from numpy.polynomial import Polynomial
from numpy import array,sqrt,real,imag
from export import ExportAllFacetings

config.useSymbolic = True

vertices.InitializeSymbolic()

def GetRoots(poly, var): #Obtains all positive real roots of a polynomial
    if poly.is_constant():
        return []
    
    Npoly = sp.N(poly)
    
    coeffs = [float(i) for i in reversed(sp.Poly(Npoly).all_coeffs())]
    f = Polynomial(coeffs)
    
    allRoots = [float(root.real) for root in f.roots() if abs(root.imag) < 10**-config.numericPrecision]
    
    return [(poly, round(root,config.numericPrecision)) for root in allRoots if root > 0]

#Parameters for the orbit types.
a = config.a
b = config.b

#Orbit types with 1 degree of freedom.
tT = vertices.GenerateVStar332(0, a, 1)
rT = vertices.GenerateVStar332(a, 1, 0)
rP = vertices.GenerateV3Star2 (0, a, 1)
tO = vertices.GenerateVStar432(0, a, 1)
tC = vertices.GenerateVStar432(a, 1, 0)
rC = vertices.GenerateVStar432(a, 0, 1)
tI = vertices.GenerateVStar532(0, a, 1)
tD = vertices.GenerateVStar532(a, 0, 1)
rD = vertices.GenerateVStar532(a, 1, 0)

#Orbit types with 2 degrees of freedom.
sT = vertices.GenerateV332(a, b, 1)
gT = vertices.GenerateVStar332(a, b, 1)
gP = vertices.GenerateV3Star2(a, b, 1)
sC = vertices.GenerateV432(a, b, 1)
gC = vertices.GenerateVStar432(a, b, 1)
sD = vertices.GenerateV532(a, b, 1)
gD = vertices.GenerateVStar532(a, b, 1)

#This expression is the determinant of
#[ p[0] p[1] p[2] 1 ]
#[ q[0] q[1] q[2] 1 ]
#[ r[0] r[1] r[2] 1 ]
#[ s[0] s[1] s[2] 1 ]
def ConfigurationEntry(p,q,r,s):
    expr = (-p[2]*q[1]*r[0] + p[1]*q[2]*r[0] + p[2]*q[0]*r[1] - p[0]*q[2]*r[1] -
             p[1]*q[0]*r[2] + p[0]*q[1]*r[2] + p[2]*q[1]*s[0] - p[1]*q[2]*s[0] - 
             p[2]*q[0]*s[1] + p[0]*q[2]*s[1] + p[1]*q[0]*s[2] - p[0]*q[1]*s[2] + 
             p[2]*r[0]*s[1] - p[1]*r[0]*s[2] - p[2]*r[1]*s[0] + p[1]*r[2]*s[0] - 
             p[0]*r[2]*s[1] + p[0]*r[1]*s[2] - q[2]*r[0]*s[1] + q[1]*r[0]*s[2] + 
             q[2]*r[1]*s[0] - q[1]*r[2]*s[0] + q[0]*r[2]*s[1] - q[0]*r[1]*s[2])
    return sp.expand(expr)

def MergePlanes(planes,newPlane):
    mergedPlanes = []
    mergingPlane = set(newPlane)
    for p in planes:
        intersectionSize = len(set(p) & set(newPlane))
        if intersectionSize >= 2:
            mergingPlane = mergingPlane.union(set(p))
        else:
            mergedPlanes.append(p)
    
    return mergedPlanes + [tuple(mergingPlane)]

def VolumeConfiguration(orbitType):
    Conf = []
    for i in range(1,len(orbitType)):
        print(i)
        for j in range(i+1,len(orbitType)):
            for k in range(j+1,len(orbitType)):
                entry = ConfigurationEntry(orbitType[0],orbitType[i],orbitType[j],orbitType[k])
                Conf.append(entry)
    return Conf

def Get1DOrbitRoots(orbitType, var):
    Conf = VolumeConfiguration(orbitType)
    roots = []
    numericalRoots = []
    
    for i in range(len(Conf)):
        if i % 100 == 0:
            print(i)
        entry = Conf[i]
        
        entryRoots = GetRoots(entry, var)
        for root in entryRoots:
            if not (root[1] in numericalRoots):
                numericalRoots.append(root[1])
                roots.append(root)
    return roots

def DeepEval1D(orbitType,var,num):
    return [array([p[0].evalf(subs={var:num}),
                   p[1].evalf(subs={var:num}),
                   p[2].evalf(subs={var:num})]) for p in orbitType]
    
def Enumerate1DOrbit(orbitType, group, var, directory, name):
    roots = Get1DOrbitRoots(orbitType, var)
    config.useSymbolic = False
    vertices.InitializeNumeric()
    
    orbits = [DeepEval1D(orbitType, var, root[1]) for root in roots]
    
    for i in range(len(orbits)):
        ExportAllFacetings(orbits[i], group, directory, name+".root"+str(roots[i][1]), ignoreTriangles=True)