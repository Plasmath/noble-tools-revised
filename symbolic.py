#Code that constructs the symbolic forms of the orbit types.
#This is only needed for the orbit types with one or more degrees of freedom
#in order to obtain exact minimal polynomials for their parameters.
import vertices
import config
import sympy as sp
import groups
from numpy.polynomial import Polynomial
from numpy import array,sqrt,real,imag
from export import ExportAllFacetings
from faceting import FindFacetings

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

def VolumeConfiguration(orbitType):
    Conf = dict()
    for i in range(1,len(orbitType)):
        for j in range(i+1,len(orbitType)):
            for k in range(j+1,len(orbitType)):
                entry = ConfigurationEntry(orbitType[0],orbitType[i],orbitType[j],orbitType[k])
                Conf[entry] = Conf.setdefault(entry,[])+[(i,j,k)]
    return Conf

def MergePlanes(planes):
    plns = [set(t) for t in sorted(planes)]
    
    mergedPlanes = []
    for p in plns:
        mergedPlanes.append(p)
        
        planesToDelete = set()
        for i in range(len(mergedPlanes)):
            for j in range(i+1,len(mergedPlanes)):
                if len(mergedPlanes[i] & mergedPlanes[j]) > 1:
                    mergedPlanes[i] = mergedPlanes[i].union(mergedPlanes[j])
                    planesToDelete.add(j)
        
        for j in reversed(list(planesToDelete)):
            del mergedPlanes[j]
    
    return mergedPlanes

#Obtain the set S_f for a given polynomial f as described in Theorem 3.21.
#This is a recursive function, and we make a minor efficiency improvement
#where if we calculate the greatest common divisor of f with a polynomial
#it is not necessary to calculate that same divisor with one of the recursive
#factors of f.
def GetCoprimeSet(factorDict, f, planes, extension, startingIndex = 0):
    
    for i in range(startingIndex,len(factorDict)):
        g = list(factorDict.keys())[i]
        res = sp.resultant(f,g)
        
        if res == 0 and f != g: #Shared common factor
            gcd = sp.gcd(f,g,extension=extension)
            gcdPlanes = MergePlanes(planes+factorDict[g])
            
            if f == gcd: #f is a factor of g
                continue
            
            if gcd.is_constant(): #The greatest common divisor was not calculated under a large enough field extension, so something went wrong
                raise Exception("Error: constant GCD encountered. The following polynomials were involved: "+str(f)+", "+str(g))
            
            h = sp.quo(f,gcd) #Computes polynomial quotient between f and gcd: the polynomial h such that f = gcd*h + r when using Euclidean division.
            #As we know already that gcd is a factor of f, r must be zero and this is the same as f/gcd.
            
            coprimeSetGCD = GetCoprimeSet(factorDict, gcd, gcdPlanes, extension, startingIndex=i)
            coprimeSetH   = GetCoprimeSet(factorDict,   h,    planes, extension, startingIndex=i)
            
            return coprimeSetGCD + coprimeSetH #Disjoint union S_gcd(f,g) U S_h
    
    return [(f,planes)]

def Copr(Conf, numVerts, extension = []):
    #Initial factoring attempt
    initialFactorsDict = dict()
    
    for c in Conf.keys():
        
        factors = [f[0] for f in sp.factor_list(c)[1]]
        
        for f in factors:
            initialFactorsDict[f] = initialFactorsDict.setdefault(f,set()).union(set(Conf[c]))
    
    keys = list(initialFactorsDict.keys())
    for c in keys:
        #Delete factors which have no negative coefficients, as these cannot have positive roots
        if all(i >= 0 for i in sp.Poly(c,a).all_coeffs()):
            del initialFactorsDict[c]
        else:
            initialFactorsDict[c] = MergePlanes(initialFactorsDict[c])
    
    initialFactors = list(initialFactorsDict.keys())
    
    filteredSets = sum((GetCoprimeSet(initialFactorsDict, f, initialFactorsDict[f], extension) for f in initialFactors),[])
    
    #The sets given may have multiple instances of the same polynomial, so we merge the planes of these together.
    coprimeFactorsDict = dict()
    
    for i in range(len(filteredSets)):
        s = filteredSets[i][0]
        newPlanes = filteredSets[i][1]
        
        #Again, we do not need to consider factors which have no negative coefficients.
        if all(i >= 0 for i in sp.Poly(s,a).all_coeffs()):
            continue
        
        currentDictValue = coprimeFactorsDict.setdefault(s, [])
        coprimeFactorsDict[s] = MergePlanes(newPlanes+currentDictValue)
    
    #Final check to make sure that no pair of polynomials shares a common factor.
    for i in coprimeFactorsDict.keys():
        for j in coprimeFactorsDict.keys():
            if sp.resultant(i,j) == 0 and i != j:
                raise Exception("Copr still has shared common factor between "+str(i)+" "+str(j))
    
    return list(coprimeFactorsDict.items())

def GetCopr(orbitType, extension = [sqrt(2)]):
    Conf = VolumeConfiguration(orbitType)
    return Copr(Conf, len(orbitType), extension = extension)

def Get1DOrbitCandidates(orbitType, copr, group):
    totalFacetings = []
    
    for s in copr:
        facetings = []
        planes = s[1]
        for p in planes:
            facetings += FindFacetings(len(orbitType), group, [0]+list(p), minCycleLength=4)
        
        if len(facetings) > 0:
            totalFacetings.append((s[0],facetings))
    
    return totalFacetings

tOcopr = GetCopr(tO, extension = [sp.sqrt(2)])
for group in groups.tOGroups:
    cand = Get1DOrbitCandidates(tO, tOcopr, group)
    print(cand)

"""
def Get1DOrbitRoots(Conf, var):
    roots = []
    numericalRoots = []
    
    for i in range(len(Conf)):
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
    
def Enumerate1DOrbit(roots, orbitType, group, var, directory, name):
    config.useSymbolic = False
    vertices.InitializeNumeric()
    
    orbits = [DeepEval1D(orbitType, var, root[1]) for root in roots]
    
    for i in range(len(orbits)):
        ExportAllFacetings(orbits[i], group, directory, name+".root"+str(roots[i][1]), ignoreTriangles=True, poly=roots[i][0])
"""