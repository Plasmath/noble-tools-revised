from symbolic import *
from faceting import FindFacetings, Generate
import sympy as sp
import groups
from export import ExportToOFF, WriteSummary

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

def GetCopr(orbitType, extension = [sp.sqrt(2)]):
    Conf = VolumeConfiguration(orbitType)
    return Copr(Conf, len(orbitType), extension = extension)

def EquivalentFaces(face, group): #Finds the set of faces containing the 0 (initial) vertex and equivalent to the given face.
    faces = [f for f in Generate(face,group) if 0 in f]
    permutedFaces = [] #We want the first vertex of our face to be the initial vertex
    for f in faces:
        i = f.index(0)
        g = f[i+1:]+f[:i]
        
        permutedFaces.append([0]+list(g))
        permutedFaces.append([0]+list(reversed(g)))
    return permutedFaces

def Get1DOrbitCandidates(orbitType, copr, group):
    totalFacetings = []
    
    for s in copr:
        facetings = []
        planes = s[1]
        for p in planes:
            facetings += FindFacetings(len(orbitType), group, [0]+list(p), minCycleLength=4)
        
        if len(facetings) > 0:
            totalFacetings.append((s[0],facetings))
    
    #Filter facetings equivalent under symmetry
    totalFacetingsFiltered = []
    for root in totalFacetings:
        faces = root[1]
        
        totalRootFacetings = []
        checkedFaces = []
        for f in faces:
            if f in checkedFaces:
                continue
            
            checkedFaces += EquivalentFaces(f, group)
            totalRootFacetings.append(f)
        
        totalFacetingsFiltered.append((root[0],totalRootFacetings))
    
    return totalFacetingsFiltered

def DeepEval1D(orbitType,num):
    return [[p[0].evalf(subs={a:num}),
             p[1].evalf(subs={a:num}),
             p[2].evalf(subs={a:num})] for p in orbitType]

def Export1DOrbitTypeFacetings(orbitType, typeCandidates, group, directory, name):
    for polyCandidates in typeCandidates:
        polynomial = polyCandidates[0]
        faces = polyCandidates[1]
        
        roots = [i.evalf() for i in sp.real_roots(polynomial)]
        for root in roots:
            if root <= 0:
                continue
            
            vertices = DeepEval1D(orbitType, root)
            
            for face in faces:
                print("Exporting faceting",face,"of",name,"at",root)
            
                ExportToOFF(vertices, face, group, directory, name+"."+str(root))
                WriteSummary(name+"."+str(root), polynomial)