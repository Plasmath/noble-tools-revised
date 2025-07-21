#symbolic.py
#Code that constructs the symbolic forms of the orbit types.
#This is only needed for the orbit types with one or more degrees of freedom
#in order to obtain exact minimal polynomials for their parameters.
import config
import sympy as sp
import groups #contains vertices.py
from numpy.polynomial import Polynomial

config.useSymbolic = True #Use symbolic mode
groups.InitializeConstants()

#Obtains all positive real roots of a polynomial.
def GetRoots(poly, var):
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
tT = groups.GenerateVStar332(0, a, 1)
rT = groups.GenerateVStar332(a, 1, 0)
rP = groups.GenerateV3Star2 (0, a, 1)
tO = groups.GenerateVStar432(0, a, 1)
tC = groups.GenerateVStar432(a, 1, 0)
rC = groups.GenerateVStar432(a, 0, 1)
tI = groups.GenerateVStar532(0, a, 1)
tD = groups.GenerateVStar532(a, 0, 1)
rD = groups.GenerateVStar532(a, 1, 0)

#Orbit types with 2 degrees of freedom.
sT = groups.GenerateV332(a, b, 1)
gT = groups.GenerateVStar332(a, b, 1)
gP = groups.GenerateV3Star2(a, b, 1)
sC = groups.GenerateV432(a, b, 1)
gC = groups.GenerateVStar432(a, b, 1)
sD = groups.GenerateV532(a, b, 1)
gD = groups.GenerateVStar532(a, b, 1)

#This expression is the determinant of
#[ p[0] p[1] p[2] 1 ]
#[ q[0] q[1] q[2] 1 ]
#[ r[0] r[1] r[2] 1 ]
#[ s[0] s[1] s[2] 1 ], 
def ConfigurationEntry(p,q,r,s):
    expr = (-p[2]*q[1]*r[0] + p[1]*q[2]*r[0] + p[2]*q[0]*r[1] - p[0]*q[2]*r[1] -
             p[1]*q[0]*r[2] + p[0]*q[1]*r[2] + p[2]*q[1]*s[0] - p[1]*q[2]*s[0] - 
             p[2]*q[0]*s[1] + p[0]*q[2]*s[1] + p[1]*q[0]*s[2] - p[0]*q[1]*s[2] + 
             p[2]*r[0]*s[1] - p[1]*r[0]*s[2] - p[2]*r[1]*s[0] + p[1]*r[2]*s[0] - 
             p[0]*r[2]*s[1] + p[0]*r[1]*s[2] - q[2]*r[0]*s[1] + q[1]*r[0]*s[2] + 
             q[2]*r[1]*s[0] - q[1]*r[2]*s[0] + q[0]*r[2]*s[1] - q[0]*r[1]*s[2])
    return sp.expand(expr)

#Obtain the volume configuration of an orbit type, with identical polynomials grouped together.
def VolumeConfiguration(orbitType):
    Conf = dict()
    for i in range(1,len(orbitType)):
        for j in range(i+1,len(orbitType)):
            for k in range(j+1,len(orbitType)):
                entry = ConfigurationEntry(orbitType[0],orbitType[i],orbitType[j],orbitType[k])
                Conf[entry] = Conf.setdefault(entry,[])+[(i,j,k)]
    return Conf

#Export the volume configuration of a given orbit type in a format easier to process with Mathematica.
def ExportConf(Conf, name, dim = 1):
    f = open("data/conf"+str(dim)+"d.txt","a")
    f.write(name+" : "+repr(Conf))
    f.write("\n")
    f.close()
    
    g = open("data/"+name+"ConfPoly.txt","w")
    for poly in Conf.keys():
        polyStr = str(poly).replace("**","^").replace("sqrt(","Sqrt[").replace(")","]").replace(" ","")
        g.write(polyStr)
        g.write("\n")
    g.close()

#Given a set of planes of coplanar vertices, group together planes which must be coplanar
#with each other. Each plane is assumed to have the vertex with index 0 as noble polyhedra
#are vertex-transitive.
def MergePlanes(planes):
    plns = [set(t) for t in sorted(tuple(i) for i in planes)]
    
    mergedPlanes = []
    for p in plns:
        mergedPlanes.append(p)
        planesToDelete = set()
        for i in range(len(mergedPlanes)):
            for j in range(i+1,len(mergedPlanes)):
                if len(mergedPlanes[i] & mergedPlanes[j]) > 1:
                    mergedPlanes[i] = mergedPlanes[i].union(mergedPlanes[j])
                    planesToDelete.add(j)
        
        for j in reversed(sorted(planesToDelete)):
            del mergedPlanes[j]
    
    return mergedPlanes

#Performs the MergePlanes operation on all planes of a dictionary.
def MergeAllPlanes(preCopr):
    Copr = preCopr.copy()
    for factor in Copr:
        tris = Copr[factor]
        Copr[factor] = MergePlanes(tris)
    return Copr

#Determine if all coefficients of a polynomial are positive.
def PositiveCoeffs(poly):
    return all(all(i >= 0 for i in sp.Poly(f,a).all_coeffs()) for f in sp.Poly(poly,b).all_coeffs())

#Given corresponding data files, imports the coprime polynomials of an
#orbit type with corresponding faceting data.
def ImportCoprData(name, Conf):
    factorsFile = open("data/"+name+"ConfFactors.txt")
    indicesFile = open("data/"+name+"ConfFactorIndices.txt")
    
    factors = []
    for i in factorsFile.readlines():
        line = i.replace("^","**").replace("Sqrt[","sp.sqrt(").replace("]",")").replace("\n","")
        factors.append(eval(line))
    
    factorPairs = [] #Each pair is (index of polynomial within Conf, index of factor within factors)
    for line in indicesFile.readlines():
        pair = tuple(int(i) for i in line.split("\t"))
        factorPairs.append(pair)
    
    factorDict = dict()
    polyList = list(Conf.keys())
    for i,j in factorPairs:
        poly = polyList[i]
        factor = factors[j]
        
        if PositiveCoeffs(factor):
            continue
        
        factorDict[factor] = factorDict.setdefault(factor, []) + Conf[poly]
    
    return factorDict