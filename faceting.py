#faceting.py
#Code that handles the process of faceting orbits and equivalence classes of orbits.

import config
import groups
from numpy import cross

#Find orbit of S under group.
def Generate(S,group):
    return [groups.Compose(g,S) for g in group]

#Obtain the set of edges of a face as a list of tuples.
def GetEdgesOfFace(face):
    return tuple(sorted([tuple(sorted([face[i],face[(i+1)%len(face)]])) for i in range(len(face))]))

def FilterPlanes(planes, group):
    filteredplanes = []
    allplanes = []
    for p in planes:
        if not (p in allplanes):
            allplanes += [tuple(sorted(q)) for q in Generate(p, group) if (0 in q)]
            filteredplanes.append(p)
    return filteredplanes

#Checks if a polyhedron is a compound given its set of faces.
def IsCompound(faces):
    numFaces = len(faces)
    
    #Construct a graph where each vertex corresponds to a face,
    #and two vertices are connected if their faces are adjacent.
    adj = {i:{j for j in range(numFaces) if len(set(faces[i]).intersection(set(faces[j]))) == 2} for i in range(numFaces)}
    
    
    #The polyhedron is a compound if and only if this graph is disconnected.
    connectedComponent = {0}
    for i in range(numFaces):
        newConnectedComponent = connectedComponent.copy()
        for j in connectedComponent:
            newConnectedComponent = newConnectedComponent.union(adj[j])
        if newConnectedComponent == connectedComponent:
            break
        connectedComponent = newConnectedComponent
    return len(connectedComponent) != numFaces

#Filters out invalid facetings which do not have an abstract polyhedral structure.
def IsValid(face,group,orbitSize):
    faces = Generate(face,group)
    
    if IsCompound(faces): #Compounds are not abstract polyhedra.
        return False
    
    edges = list(GetEdgesOfFace(face))
    combinedEdges = list({GetEdgesOfFace(face) for face in faces})
    totalEdges = sum((list(i) for i in combinedEdges), [])
    
    #All edges must be adjacent to exactly 2 faces.
    for e in edges:
        if totalEdges.count(e) != 2: #checks if this edge appears in exactly 2 faces
            return False
    
    return True

#Given a group and plane of vertices, finds all polygons in that plane which generate
#a valid noble polyhedron. minCycleLength gives the minimum number of sides of a polygon,
#which can be increased in certain situations for a performance boost.
def FindFacetings(orbitSize, group, plane, minCycleLength = 3): #Find facetings of an orbit within a given plane.
    #NOTE: This may find certain extra invalid facetings which do not have two faces to an edge.
    #This is due to some edges in the plane requiring other edges in order to exist, which is not implemented.  

    planes = Generate(plane,group)
    
    graph = set([]) #Adjacency graph.
    for p in planes:
        intersection = tuple(set(p) & set(plane))
        if len(intersection) == 2:
            graph.add(intersection)
    
    #We store the graph in a different manner so we may look at a given vertex
    #when finding cycles within the graph. In particular,
    #adj[i] gives the indices of the vertices adjacent to i
    adj = {i:{j for j in plane if (i,j) in graph or (j,i) in graph} for i in plane}
    
    #the cycles must contain the main vertex of the graph, so at least 2
    #edges must exit this vertex for there to be facetings
    if len(adj[0]) < 2:
        return []
    
    #perform BFS to find cycles containing the 0 point within the graph
    partialCycles = [[0]]
    cycles = []
    for i in range(len(plane)): #maximum cycle length contains all points in graph
        newPartialCycles = []
        
        #for each partial cycle, add all possible next vertices in the cycle
        #and add it to the next iteration of partialCycles.
        for partial in partialCycles:
            possibleNextVertices = adj[partial[-1]]
            
            for v in possibleNextVertices:
                if v == 0 and len(partial) >= minCycleLength: #the cycle has been completed, so add to cycles
                    cycles.append(partial)
                elif not (v in partial): #check that cycle does not wrap around on itself
                    newPartialCycles.append(partial + [v])
        
        partialCycles = newPartialCycles.copy() #prepare for next iteration
    
    #filter out cycles which are the same but orientation-reversed
    filtered = [cycle for cycle in cycles if cycle[1] < cycle[-1]]
    
    #return only valid facetings
    return [cycle for cycle in filtered if IsValid(cycle, group, orbitSize)]

#Given an orbit and three points from the orbit, find all points of the orbit which are
#in the plane defined by the three points.
def GetPlane(ind1,ind2,ind3,orbit):
    p1,p2,p3 = orbit[ind1],orbit[ind2],orbit[ind3]
    prec = 10**(-config.numericPrecision+4)
    c = cross(p2 - p1, p3 - p1)
    return [i for i in range(len(orbit)) if abs((orbit[i]-p1).dot(c)) < prec]

#Given an orbit, find all facetings that exist in any plane (as FindFacetings only does this in a single plane).
def FacetAll(orbit, group, ignoreTriangles = False): #Find all facetings of an orbit in all possible planes
    orbitSize = len(orbit)    

    if ignoreTriangles:
        minCycleLength = 4
    else:
        minCycleLength = 3

    #obtain set of all planes
    planes = set([])
    for i in range(1,len(orbit)):
        for j in range(i+1,len(orbit)):
            planes.add(tuple(GetPlane(0,i,j,orbit)))
    
    #find facetings in each plane
    facetings = []
    for plane in FilterPlanes(planes, group):
        facetings += FindFacetings(orbitSize, group, plane, minCycleLength = minCycleLength)
    return facetings

#In an orbit type with more than one degree of freedom, find the facetings (if any)
#of the minimum equivalence class.
def FacetMinimalEquivalenceClass(orbitType, sharedPlanes, group):
    totalFacetings = []
    orbitSize = len(orbitType)
    
    #Check for noble polyhedra with triangular faces
    for i in range(1,orbitSize):
        for j in range(i+1,orbitSize):
            cycle = [0,i,j]
            if IsValid(cycle, group, orbitSize):
                totalFacetings.append(cycle)
    
    #Check for noble polyhedra in shared planes - triangles have already been checked, so we let minCycleLength = 4.
    for plane in sharedPlanes:
        totalFacetings += FindFacetings(orbitSize, group, [0]+list(plane), minCycleLength=4)
    
    return totalFacetings