import groups
import vertices
import config
from numpy import cross

def Generate(S,group):
    return [groups.Compose(g,S) for g in group]

def PlaneStabilizer(plane,group):
    for g in group:
        print(list(sorted(groups.Compose(g,plane))), list(sorted(plane)))
    
    return [g for g in group if list(sorted(groups.Compose(g,plane)))==list(sorted(plane))]

def IsValid(face,stabilizer,plane):
    edges = [list(sorted([face[i],face[(i+1)%len(face)]])) for i in range(len(face))]
    #print(edges)
    
    for e in edges:
        iterator = [(list(sorted(i)) in edges) for i in Generate(e,stabilizer)]
        #print(iterator, any(iterator), all(iterator))
        if any(iterator) and not all(iterator):
            return False
    return True

def FindFacetings(orbit, group, plane): #Find facetings of an orbit within a given plane.
    #NOTE: This may find certain nonpolyhedra which contain coplanar faces sharing an edge.
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
                if v == 0 and len(partial) > 2: #the cycle has been completed, so add to cycles
                    cycles.append(partial)
                elif not (v in partial): #check that cycle does not wrap around on itself
                    newPartialCycles.append(partial + [v])
        
        partialCycles = newPartialCycles.copy() #prepare for next iteration
    
    #filter out cycles which are the same but orientation-reversed
    filtered = [cycle for cycle in cycles if cycle[1] < cycle[-1]]
    
    #return only valid facetings
    stabilizer = PlaneStabilizer(plane, group)
    return [cycle for cycle in filtered if IsValid(cycle, stabilizer, plane)]

def GetPlane(ind1,ind2,ind3,orbit):
    p1,p2,p3 = orbit[ind1],orbit[ind2],orbit[ind3]
    prec = 10**-config.numericPrecision
    c = cross(p2 - p1, p3 - p1)
    return [i for i in range(len(orbit)) if abs((orbit[i]-p1).dot(c)) < prec]

def FacetAll(orbit, group): #Find all facetings of an orbit in all possible planes
    #obtain set of all planes
    planes = set([])
    for i in range(1,len(orbit)):
        for j in range(i+1,len(orbit)):
            planes.add(tuple(GetPlane(0,i,j,orbit)))
    
    #find facetings in each plane
    facetings = []
    for plane in planes:
        facetings += FindFacetings(orbit, group, plane)
    return facetings
    
V = vertices.GenerateVStar532(1, 0, 0)
f = FindFacetings(V, groups.DGroupStar532, [0, 3, 13, 4, 19, 14])

print([tuple(V[i]) for i in [0, 3, 4, 13, 14, 15]])


#print(f,len(f))