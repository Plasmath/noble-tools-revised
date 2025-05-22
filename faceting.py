import groups
import config
from numpy import cross

def Generate(S,group):
    return [groups.Compose(g,S) for g in group]

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

def IsCompound(faces): #Checks if a polyhedron is a compound
    numFaces = len(faces)
    
    adj = {i:{j for j in range(numFaces) if len(set(faces[i]).intersection(set(faces[j]))) == 2} for i in range(numFaces)}
    
    connectedComponent = {0}
    for i in range(numFaces):
        newConnectedComponent = connectedComponent.copy()
        for j in connectedComponent:
            newConnectedComponent = newConnectedComponent.union(adj[j])
        if newConnectedComponent == connectedComponent:
            break
        connectedComponent = newConnectedComponent
    return len(connectedComponent) != numFaces

def IsValid(face,group,orbitSize): #Removes invalid facetings which do not have 2 faces to an edge.
    faces = Generate(face,group)
    
    if IsCompound(faces):
        return False
    
    edges = list(GetEdgesOfFace(face))
    combinedEdges = list({GetEdgesOfFace(face) for face in faces})
    totalEdges = sum((list(i) for i in combinedEdges), [])
    
    for e in edges:
        if totalEdges.count(e) != 2: #checks if this edge appears in exactly 2 faces
            return False
    
    return True

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

def GetPlane(ind1,ind2,ind3,orbit):
    p1,p2,p3 = orbit[ind1],orbit[ind2],orbit[ind3]
    prec = 10**(-config.numericPrecision+4)
    c = cross(p2 - p1, p3 - p1)
    return [i for i in range(len(orbit)) if abs((orbit[i]-p1).dot(c)) < prec]

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
