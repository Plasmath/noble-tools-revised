import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from faceting import Generate

def GeneratePolyhedronPlot(vertices,edges,face,name):
    faceCoordinates = np.array([ [vertices[i] for i in face+[face[0]]] ])
    edgeCoordinates = np.array([ [vertices[e[0]],vertices[e[1]]] for e in edges ])
    
    ax = plt.figure().add_subplot(projection='3d')
    
    edgeCollection = Line3DCollection(edgeCoordinates, linewidth=0.5)
    
    ax.add_collection3d(edgeCollection)
    
    faceCollection = Poly3DCollection(faceCoordinates, alpha=1, color="gray")
    ax.add_collection3d(faceCollection)
    
    
    #ax.set_axis_off()
    #ax.set_proj_type('persp', focal_length=0.1)
    
    plt.savefig('images/'+name+'.png')

def MakePlot(vertices,face,group,name):
    edges = set([])
    L = len(face)
    
    for f in Generate(face,group):
        faceEdges = set(tuple(sorted([f[i],f[(i+1)%L]])) for i in range(L))
        edges = edges.union(faceEdges)
    
    GeneratePolyhedronPlot(vertices, list(edges), face, name)