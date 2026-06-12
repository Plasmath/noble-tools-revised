#symbolic2D.py
#Code specific to the enumeration of noble polyhedra found within
#orbit types with 2 degrees of freedom.

from symbolic import *
from faceting import FindFacetings, Generate
from sympy import sqrt, resultant
from export import ExportToOFF, WriteSummary

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
def SelectiveMergePlanes(mainPlanes, sharedPlanes, includeOtherPlanes = False):
    totalPlanes = mainPlanes.copy()
    remainingPlanes = []
    for s in sharedPlanes:
        if any( len(s & i) > 1 for i in mainPlanes ):
            totalPlanes.append(s)
        elif includeOtherPlanes:
            remainingPlanes.append(s)
    return MergePlanes(totalPlanes) + remainingPlanes

#Determine the full set of planes for all the different coprime polynomials.
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
    
    pairDataDictA = dict() #For information about how this dictionary is structured, see enumerate2D-D.py.
    pairDataDictB = dict()
    
    for line in f.readlines():
        split = line.split("\t") #Extract data from line
        
        pair      = (int(split[0]),int(split[1])) #Pair of polynomials being analyzed for intersections.
        var       = split[2]                      #Variable for which a root has been found.
        root      = float(split[3].split("`")[0]) #The actual approximate value of the root.
        precision = float(split[3].split("`")[1]) #The known precision of the value given.
        
        pairDataDictA[pair] = pairDataDictA.setdefault(pair, [])
        pairDataDictB[pair] = pairDataDictB.setdefault(pair, [])
        if var == "a":
            pairDataDictA[pair].append( (root,precision) )
        else:
            pairDataDictB[pair].append( (root,precision) )
    
    #Filter situations where there are no valid values for one of a or b,
    #in which case we do not need to consider the existence of such an intersection.
    for key in list(pairDataDictA.keys()):
        if len(pairDataDictA[key]) == 0 or len(pairDataDictB[key]) == 0:
            del pairDataDictA[key]
            del pairDataDictB[key]
    
    return pairDataDictA, pairDataDictB

#Returns True if the upper and lower bounds of n1 and n2 overlap.
def AreClose(n1,n2):
    uncertainty = 10**(len(str(int(n1[0])))-int(n1[1])) + 10**(len(str(int(n2[0])))-int(n2[1]))    
    return min(n1[0],n2[0]) + uncertainty >= max(n1[0],n2[0])

#Find the connected component of a graph that contains a given node.
def GetConnectedComponent(adjDict,startNode):
    connectedComponent = [startNode] #List of points in connected component.
    
    nodeConnections = adjDict[startNode] #List of points that still need to be searched for additional graph connections
    del adjDict[startNode]
    while len(nodeConnections) > 0:
        for i in nodeConnections.copy():
            nodeConnections.remove(i)
            
            if i in adjDict:
                connectedComponent.append(i)
                nodeConnections += adjDict[i]
                del adjDict[i]
    
    return connectedComponent

def GetParameterValues(pairData):
    aValues = dict()
    bValues = dict()
    
    for pair in pairData[0].keys():
        roots = pairData[0][pair]
        for root in roots:
            aValues[root] = aValues.setdefault(root,[]) + [pair]
    for pair in pairData[1].keys():
        roots = pairData[1][pair]
        for root in roots:
            bValues[root] = bValues.setdefault(root,[]) + [pair]
    
    return aValues, bValues

def ImportIntersectionData(name):
    data = open("data/"+name+"Intersections.txt").readlines()
    data = [eval(i.replace("{","[").replace("}","]")) for i in data]
    
    intersectionDict = dict()
    for pair, intersections in data:
        for intersection in intersections:
            intersectionDict[tuple(intersection)] = intersectionDict.setdefault(tuple(intersection), []) + [pair]
    
    return intersectionDict

def ImportFacetingData(name):
    polynomialsFile = open("data/"+name+"CoprList.txt")
    coprFile = open("data/"+name+"Copr.txt").readlines()
    copr = eval(coprFile[0])
    sharedPlanes = eval(coprFile[1])
    
    polynomials = []
    for i in polynomialsFile.readlines():
        line = i.replace("^","**").replace("Sqrt[2]","sp.sqrt(2)").replace("Sqrt[5]","sp.sqrt(5)").replace("\n","")
        polynomials.append(eval(line))
    
    return polynomials, copr, sharedPlanes

def DetermineIntersectionData(intersections, poly, copr, sharedPlanes):
    intersectionData = dict([])
    keys = list(intersections.keys())
    for n in range(len(keys)):
        k = keys[n]
        
        pairs = intersections[k]
        curves = set(sum(pairs, []))
        
        mergedPlanes = []
        for i in curves:
            mergedPlanes = MergePlanes(mergedPlanes + copr[poly[i-1]])
        
        intersectionData[k] = SelectiveMergePlanes(mergedPlanes, sharedPlanes, includeOtherPlanes=True)
        
    return intersectionData

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

def Get2DOrbitTypeFacetings(orbitType, planesDictionary, group):
    totalFacetings = []
    
    for i,j in planesDictionary.keys():
        facetings = []
        planes = planesDictionary[(i,j)]
        for p in planes:
            facetings += FindFacetings(len(orbitType), group, [0]+list(p), minCycleLength=4)
        
        if len(facetings) > 0:
            totalFacetings.append((i,j,facetings))
    
    #Filter facetings equivalent under symmetry
    totalFacetingsFiltered = []
    for root in totalFacetings:
        faces = root[2]
        
        totalRootFacetings = []
        checkedFaces = []
        for f in faces:
            if f in checkedFaces:
                continue
            
            checkedFaces += EquivalentFaces(f, group)
            totalRootFacetings.append(f)
        
        totalFacetingsFiltered.append((root[0],root[1],totalRootFacetings))
    
    return totalFacetingsFiltered

def ImportAtlasData(name):
    atlasFile = open("data/"+name+"Atlas.txt")
    
    lines = atlasFile.readlines()
    
    aValues = [float(i.split("`")[0]) for i in lines[0][1:-1].split(",")]
    bValues = [float(i.split("`")[0]) for i in lines[3][1:-1].split(",")]
    
    aPolynomials = []
    for i in lines[1][1:-2].split(","):
        aPolynomials.append( eval(i[:-4].replace("^","**").replace("Sqrt[2]","sp.sqrt(2)").replace("Sqrt[5]","sp.sqrt(5)").replace("Sqrt[10]","sp.sqrt(10)")) )
    
    bPolynomials = []
    for i in lines[4][1:-2].split(","):
        bPolynomials.append( eval(i[:-4].replace("^","**").replace("Sqrt[2]","sp.sqrt(2)").replace("Sqrt[5]","sp.sqrt(5)").replace("Sqrt[10]","sp.sqrt(10)")) )
    
    atlas = [aValues, aPolynomials, bValues, bPolynomials]
    
    return atlas

#Substitutes the values of two constants for the parameters of
#a given orbit type with 2 degree of freedom.
def DeepEval2D(orbitType,aVal,bVal):
    return [[p[0].evalf(subs={a:aVal, b:bVal}),
             p[1].evalf(subs={a:aVal, b:bVal}),
             p[2].evalf(subs={a:aVal, b:bVal})] for p in orbitType]

def Export2DOrbitTypeFacetings(orbitType, atlas, facetings, group, directory, name):
    aValues, aPolynomials, bValues, bPolynomials = atlas
    
    for vertexSet in facetings:
        aAtlasIndex, bAtlasIndex, faces = vertexSet
        
        aValue = aValues[aAtlasIndex - 1]
        aPoly = aPolynomials[aAtlasIndex - 1]
        bValue = bValues[bAtlasIndex - 1]
        bPoly = bPolynomials[bAtlasIndex - 1]
        
        vertices = DeepEval2D(orbitType, aValue, bValue)
        
        for i in range(len(faces)):
            face = faces[i]
            print("Exporting faceting",face,"of",name,"at",[aValue,bValue])
            
            fileName = name+"."+str(aValue)+"."+str(bValue)+"."+str(i)
            
            ExportToOFF(vertices, face, group, directory, fileName)
            WriteSummary(fileName, [aPoly, bPoly], secondParameter=True)
        