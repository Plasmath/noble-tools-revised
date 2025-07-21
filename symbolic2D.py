#symbolic2D.py
#Code specific to the enumeration of noble polyhedra found within
#orbit types with 2 degrees of freedom.

from symbolic import *
from faceting import FindFacetings, Generate
import sympy as sp
from export import ExportToOFF, WriteSummary
import time

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
            