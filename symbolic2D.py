#symbolic2D.py
#Code specific to the enumeration of noble polyhedra found within
#orbit types with 2 degrees of freedom.

from symbolic import *
from faceting import FindFacetings
from sympy import sqrt, resultant

#Given the appropriate data files, imports the coprime polynomials
#for a given orbit type and collects the necessary information to
#obtain the relevant faceting data.
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

#Exports the coprime polynomials and faceting data for a given orbit type.
def ExportCopr(name, copr, sharedPlanes):
    f = open("data/"+name+"Copr.txt", "w")
    f.write(repr(copr))
    f.write("\n")
    f.write(repr(sharedPlanes))
    f.close()

#Merges two sets of planes together, with the condition that the
#a plane from the second set is only added to the final list if
#it will merge with a plane from the first set.
def SelectiveMergePlanes(mainPlanes, sharedPlanes):
    totalPlanes = mainPlanes.copy()
    for s in sharedPlanes:
        if any( len(s & i) > 1 for i in mainPlanes ):
            totalPlanes.append(s)
    return MergePlanes(totalPlanes)

def MergeAllWithShared(copr, sharedPlanes):
    for s in copr.keys():
        copr[s] = SelectiveMergePlanes(copr[s], sharedPlanes)
    return copr

#Test a given set of planes for abstract noble facetings under
#a given symmetry group.
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

#Determines if two sets of planes contain two planes that would
#merge to form a larger plane.
def RequiresMerge(sPlanes,tPlanes):
    for p1 in sPlanes:
        if any((len(p2 & p1) > 1 and p2 != p1) for p2 in tPlanes):
            return True
    return False

#Finds all pairs of coprime polynomials within an orbit type
#such that the critical planes of these polynomials will merge
#to form larger planes at their intersections.
def FindCriticalPairs(copr):
    factors = list(copr.keys())
    criticalPairs = []
    for i in range(len(factors)):
        sPlanes = copr[factors[i]]
        for j in range(i+1,len(factors)):
            tPlanes = copr[factors[j]]
            if RequiresMerge(sPlanes,tPlanes):
                criticalPairs.append((i,j))
    return criticalPairs

def ExportCriticalPairs(name, pairs, coprList):
    f = open("data/"+name+"CoprList.txt", "w")
    g = open("data/"+name+"Pairs.txt", "w")
    
    for poly in coprList: #Reformat polynomials into Mathematica's format.
        f.write(str(poly).replace("**","^").replace("sqrt(2)","Sqrt[2]").replace("sqrt(5)","Sqrt[5]").replace(" ",""))
        f.write("\n")
    f.close()
    
    for pair in pairs:
        g.write(str(pair[0]+1)) #Add 1 because Mathematica is a 1-indexing language.
        g.write(" ")
        g.write(str(pair[1]+1))
        g.write("\n")
    g.close()

def ImportPairData(name):
    f = open("data/"+name+"RootData.txt")
    
    pairDataDict = dict() #For information about how this dictionary is structured, see enumerate2D-D.py.
    
    for line in f.readlines():
        split = line.split("\t") #Extract data from line
        
        pair      = (int(split[0]),int(split[1])) #Pair of polynomials being analyzed for intersections.
        var       = split[2]                      #Variable for which a root has been found.
        root      = float(split[3].split("`")[0]) #The actual approximate value of the root.
        precision = float(split[3].split("`")[1]) #The known precision of the value given.
        
        pairDataDict[pair] = pairDataDict.setdefault(pair, [[],[]])
        if var == "a":
            pairDataDict[pair][0].append( (root,precision) )
        else:
            pairDataDict[pair][1].append( (root,precision) )
    
    #Filter situations where there are no valid values for one of a or b,
    #in which case we do not need to consider the existence of such an intersection.
    for key in list(pairDataDict.keys()):
        data = pairDataDict[key]
        if len(data[0]) == 0 or len(data[1]) == 0:
            del pairDataDict[key]
    
    return pairDataDict

def AreClose(n1,n2):
    uncertainty = 10**(len(str(int(n1[0])))-int(n1[1])) + 10**(len(str(int(n2[0])))-int(n2[1]))    
    return min(n1[0],n2[0]) + uncertainty >= max(n1[0],n2[0])
