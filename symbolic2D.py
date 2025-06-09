from symbolic import *
from faceting import FindFacetings, Generate
import sympy as sp
from export import ExportToOFF, WriteSummary
import time

def PositiveCoeffs(poly):
    return all(all(i >= 0 for i in sp.Poly(f,a).all_coeffs()) for f in sp.Poly(poly,b).all_coeffs())

def Copr(Conf, numVerts, extension = []):
    startTime = time.time()
    
    #Initial factoring attempt
    initialFactorsDict = dict()
    
    sharedPlanes = [tuple(i) for i in MergePlanes(Conf[0])]#Planes shared by all members of the orbit type.
    
    for c in Conf.keys():
        factors = [f[0] for f in sp.factor_list(c)[1]]
        
        for f in factors:
            #Delete factors which have no negative coefficients, as these cannot have positive roots
            if PositiveCoeffs(f):
                continue
            initialFactorsDict[f] = initialFactorsDict.setdefault(f,set()).union(set(Conf[c]))
    
    initTime = time.time()
    print("Factor initialization:",initTime-startTime)
    
    keys = list(initialFactorsDict.keys())
    for c in keys:
        initialFactorsDict[c] = MergePlanes(initialFactorsDict[c])
    
    initialFactors = list(initialFactorsDict.keys())
    
    filteredSets = []
    for i in range(len(initialFactors)):
        print(i,"/",len(initialFactors))
        f = initialFactors[i]
        filteredSets += GetCoprimeSet(initialFactorsDict, f, initialFactorsDict[f], extension)
    
    filterTime = time.time()
    print("Filtering:",filterTime-initTime)
    
    #The sets given may have multiple instances of the same polynomial, so we merge the planes of these together.
    coprimeFactorsDict = dict()
    
    for i in range(len(filteredSets)):
        s = filteredSets[i][0]
        newPlanes = filteredSets[i][1]
        
        #Again, we do not need to consider factors which have no negative coefficients.
        if PositiveCoeffs(s):
            continue
        
        currentDictValue = coprimeFactorsDict.setdefault(s, sharedPlanes)
        coprimeFactorsDict[s] = MergePlanes(newPlanes+currentDictValue)
    
    #Final check to make sure that no pair of polynomials shares a common factor.
    for i in coprimeFactorsDict.keys():
        for j in coprimeFactorsDict.keys():
            if sp.resultant(i,j) == 0 and i != j:
                raise Exception("Copr still has shared common factor between "+str(i)+" "+str(j))
    
    print("Final check:",time.time()-filterTime)
    
    return list(coprimeFactorsDict.items())

def GetCopr(orbitType, Conf, extension = [sp.sqrt(2)]):
    return Copr(Conf, len(orbitType), extension = extension)

def ExportCopr(Copr, name):
    f = open("data/copr.txt","a")
    f.write(name+" : "+repr(Copr))
    f.write("\n")
    f.close()