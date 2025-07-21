#symbolic1D.py
#Code specific to the enumeration of orbit types with 1 degree of freedom. 

from sympy import sqrt
from symbolic import *
from faceting import FindFacetings, Generate
from export import ExportToOFF, WriteSummary

#Obtain coprime polynomials and planes of an orbit type with 1 degree of freedom,
#given the volume configruation. In all cases factorization under a field extension
#is required, which is included here.
def Copr(Conf, numVerts, extension = []):
    #Initial factoring attempt
    initialFactorsDict = dict()
    
    sharedPlanes = set() #Planes shared by all members of the orbit type.
    
    for c in Conf.keys():
        if c == 0:
            sharedPlanes = sharedPlanes.union(set(Conf[c]))
        
        factors = [f[0] for f in sp.factor_list(c)[1]]
        
        for f in factors:
            initialFactorsDict[f] = initialFactorsDict.setdefault(f,set()).union(set(Conf[c]))
    
    sharedPlanes = MergePlanes(sharedPlanes)
    
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
        
        currentDictValue = coprimeFactorsDict.setdefault(s, sharedPlanes)
        coprimeFactorsDict[s] = MergePlanes(newPlanes+currentDictValue)
    
    #Final check to make sure that no pair of polynomials shares a common factor.
    for i in coprimeFactorsDict.keys():
        for j in coprimeFactorsDict.keys():
            if sp.resultant(i,j) == 0 and i != j:
                raise Exception("Copr still has shared common factor between "+str(i)+" "+str(j))
    
    return list(coprimeFactorsDict.items())

#Slightly more streamlined and convenient version of the Copr function.
def GetCopr(orbitType, extension = [sp.sqrt(2)]):
    Conf = VolumeConfiguration(orbitType)
    return Copr(Conf, len(orbitType), extension = extension)

#Finds the set of faces containing the 0 (initial) vertex and equivalent to the given face.
def EquivalentFaces(face, group):
    faces = [f for f in Generate(face,group) if 0 in f]
    permutedFaces = [] #We want the first vertex of our face to be the initial vertex
    for f in faces:
        i = f.index(0)
        g = f[i+1:]+f[:i]
        
        permutedFaces.append([0]+list(g))
        permutedFaces.append([0]+list(reversed(g)))
    return permutedFaces

#Abstractly facets coprime polynomials to check for possible noble facetings.
#From there, it may be checked whether such a faceting can actually be realized.
def Get1DOrbitCandidates(orbitType, copr, group):
    totalFacetings = []
    
    for s in copr:
        facetings = []
        planes = copr[s]
        for p in planes:
            facetings += FindFacetings(len(orbitType), group, [0]+list(p), minCycleLength=4)
        
        if len(facetings) > 0:
            totalFacetings.append((s,facetings))
    
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

#Substitutes the values of a constant for the parameter of
#a given orbit type with 1 degree of freedom.
def DeepEval1D(orbitType,num):
    return [[p[0].evalf(subs={a:num}),
             p[1].evalf(subs={a:num}),
             p[2].evalf(subs={a:num})] for p in orbitType]

#Given an orbit type and output from Get1DOrbitCandidates,
#determines the possible noble realizations of these abstract
#polyhedra and exports them as .off files.
def Export1DOrbitTypeFacetings(orbitType, typeCandidates, group, directory, name):
    for polyCandidates in typeCandidates:
        polynomial = polyCandidates[0]
        faces = polyCandidates[1]
        
        roots = [i.evalf() for i in sp.real_roots(polynomial.evalf())]
        for root in roots:
            if root <= 0:
                continue
            
            vertices = DeepEval1D(orbitType, root)
            
            for i in range(len(faces)):
                face = faces[i]
                print("Exporting faceting",face,"of",name,"at",root)
            
                ExportToOFF(vertices, face, group, directory, name+"."+str(root)+"."+str(i))
                WriteSummary(name+"."+str(root), polynomial)

