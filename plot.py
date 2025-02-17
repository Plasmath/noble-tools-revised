import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from faceting import Generate
from matplotlib import cm
from matplotlib.colors import LightSource

def GeneratePolyhedronPlot(vertices,edges,face,name):
    faceCoordinates = np.array([ [vertices[i] for i in face+[face[0]]] ])
    edgeCoordinates = np.array([ [vertices[e[0]],vertices[e[1]]] for e in edges ])
    
    norms = np.array([e[0][0]*e[0][0]+e[0][1]*e[0][1]+e[0][2]*e[0][2]+
             e[1][0]*e[1][0]+e[1][1]*e[1][1]+e[1][2]*e[1][2] for e in edgeCoordinates])
    ls = LightSource(270, 45)
    rgb = ls.shade(norms, cmap=cm.gist_earth)
    
    ax = plt.figure().add_subplot(projection='3d')
    
    edgeCollection = Line3DCollection(edgeCoordinates, cmap=cm.coolwarm, linewidth=0.5, facecolors=rgb)
    
    ax.add_collection3d(edgeCollection)
    
    faceCollection = Poly3DCollection(faceCoordinates, alpha=.7, color="gray")
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