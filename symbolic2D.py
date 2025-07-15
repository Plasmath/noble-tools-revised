from symbolic import *
from faceting import FindFacetings, Generate
import sympy as sp
from export import ExportToOFF, WriteSummary
import time

def PositiveCoeffs(poly):
    return all(all(i >= 0 for i in sp.Poly(f,a).all_coeffs()) for f in sp.Poly(poly,b).all_coeffs())

factorDict = dict()
def GetCoprimeSet(f, planes, extension, startingIndex = 0):
    global factorDict
    
    if PositiveCoeffs(f):
        return []
    
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
            
            coprimeSetGCD = GetCoprimeSet(gcd, gcdPlanes, extension, startingIndex=i)
            coprimeSetH   = GetCoprimeSet(  h,    planes, extension, startingIndex=i)
            
            return coprimeSetGCD + coprimeSetH #Disjoint union S_gcd(f,g) U S_h
    
    return [(f,planes)]

def Copr(Conf, numVerts, extension = []):
    startTime = time.time()
    
    #Initial factoring attempt
    initialFactorsDict = dict()
    
    if 0 in Conf:
        sharedPlanes = [tuple(i) for i in MergePlanes(Conf[0])]#Planes shared by all members of the orbit type.
    else:
        sharedPlanes = []
    
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
    
    global factorDict
    factorDict = initialFactorsDict
    
    filteredSets = []
    for i in range(len(initialFactors)):
        f = initialFactors[i]
        print(i,"/",len(initialFactors))
        filteredSets += GetCoprimeSet(f, initialFactorsDict[f], extension)
    
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
        
    print("Total:",time.time()-startTime)
    
    return list(coprimeFactorsDict.items())

def GetCopr(orbitType, Conf, extension = [sp.sqrt(2)]):
    return Copr(Conf, len(orbitType), extension = extension)

def ExportConf(Conf, name):
    f = open("data/conf.txt","a")
    f.write(name+" : "+repr(Conf))
    f.write("\n")
    f.close()
    
    g = open("data/"+name+"ConfPoly.txt","w")
    for poly in Conf.keys():
        polyStr = str(poly).replace("**","^").replace("sqrt(","Sqrt[").replace(")","]").replace(" ","")
        g.write(polyStr)
        g.write("\n")
    g.close()
    
def ImportCoprData(name, Conf):
    factorsFile = open("data/"+name+"ConfFactors.txt")
    indicesFile = open("data/"+name+"ConfFactorIndices.txt")
    
    factors = []
    for i in factorsFile.readlines():
        line = i.replace("^","**").replace("Sqrt[2]","sp.sqrt(2)").replace("Sqrt[5]","sp.sqrt(5)").replace("\n","")
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

def MergeAllPlanes(preCopr):
    Copr = preCopr.copy()
    for factor in Copr:
        tris = Copr[factor]
        Copr[factor] = MergePlanes(tris)
    return Copr

def ExportCopr(name, copr, sharedPlanes):
    f = open("data/"+name+"Copr.txt", "w")
    f.write(repr(copr))
    f.write("\n")
    f.write(repr(sharedPlanes))
    f.close()

def SelectiveMergePlanes(mainPlanes, sharedPlanes):
    totalPlanes = mainPlanes.copy()
    for s in sharedPlanes:
        if any( len(s & i) > 1 for i in mainPlanes ):
            totalPlanes.append(s)
    return MergePlanes(totalPlanes)

def Get2DCoprCandidates(orbitType, copr, sharedPlanes, group,):
    orbitSize = len(orbitType)
    totalFacetingCandidates = []
    
    factors = list(copr.keys())
    for i in range(len(factors)):
        s = factors[i]
        sFacetings = []
        
        planes = SelectiveMergePlanes(copr[s],sharedPlanes)
        
        for plane in planes:
            sFacetings += FindFacetings(orbitSize, group, [0]+list(plane), minCycleLength=4)
        totalFacetingCandidates += [(s,cycle) for cycle in sFacetings]
    
    return totalFacetingCandidates

ImportCopr = lambda name : [eval(i) for i in open("data/"+name+"Copr.txt").readlines()]

def RequiresMerge(sPlanes,tPlanes):
    for p1 in sPlanes:
        if any((len(p2 & p1) > 1 and p2 != p1) for p2 in tPlanes):
            return True
    return False

def FindCriticalPairs(copr, sharedPlanes):
    factors = list(copr.keys())
    criticalPairs = []
    for i in range(len(factors)):
        sPlanes = SelectiveMergePlanes(copr[factors[i]], sharedPlanes)
        for j in range(i+1,len(factors)):
            tPlanes = SelectiveMergePlanes(copr[factors[j]], sharedPlanes)
            if RequiresMerge(sPlanes,tPlanes):
                criticalPairs.append((i,j))
    return criticalPairs
            