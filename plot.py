import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from faceting import Generate
from vertices import GenerateVStar432
from groups import tOGroupStar432

tOedges = [(0,1),(0,2),(0,3),(1,5),(1,6),(2,4),(2,7),(3,6),(3,9),(4,8),(4,10),(5,8),(5,11),(6,13),(7,10),(7,12),(8,14),(9,12),(9,15),(10,17),(11,14),(11,16),(12,18),(13,16),(13,19),(14,20),(15,18),(15,19),(16,21),(17,20),(17,22),(18,22),(19,21),(20,23),(21,23),(22,23)]

def GeneratePolyhedronPlot(vertices,edges,face,name):
    faceCoordinates = np.array([ [vertices[i] for i in face+[face[0]]] ])
    edgeCoordinates = np.array([ [vertices[e[0]],vertices[e[1]]] for e in edges ])
    
    ax = plt.figure().add_subplot(projection='3d')
    
    sliceVal = 10
    
    edgeCollection1 = Line3DCollection(edgeCoordinates[10:], linewidth=0.5, linestyle="--", color="black")
    
    ax.add_collection3d(edgeCollection1)
    
    faceCollection = Poly3DCollection(faceCoordinates, alpha=1, color="gray")
    ax.add_collection3d(faceCollection)
    
    #ax.set_axis_off()
    ax.set_proj_type('persp', focal_length=0.5)
    
    plt.savefig('images/'+name+'.png')

def MakePlot(vertices,edges,face,group,name):
    GeneratePolyhedronPlot(vertices, set(edges), face, name)