from faceting import Generate, FacetAll
from plot import MakePlot
from sympy import real_roots

summary = open("3dmodels/summary.txt","w")

summary.write(
"""NOTE: These polynomials are not simplified. In some cases,
there may be a simpler factor of the given polynomial that gives the same root.

""")

summary.close()

def ExportToOFF(vertices, face, group, directory, name):
    #faces with symmetry will have duplicates under Generate.
    #This filters those extra faces out.
    faces = []
    facesets = []
    for f in Generate(face,group):
        if not (set(f) in facesets):
            faces.append(f)
            facesets.append(set(f))
    
    file = open(directory+"/"+name+".off", "w")
    
    file.write("OFF\n")
    file.write(str(len(vertices))+" "+str(len(faces))+" 0\n\n")
    
    for v in vertices:
        file.write(str(v[0])+" "+str(v[1])+" "+str(v[2])+"\n")
    
    file.write("\n")
    
    for f in faces:
        file.write(str(len(f))+" "+" ".join(str(i) for i in f)+"\n")
    
    #MakePlot(vertices, face, group, name)
        
def WriteSummary(name, polynomial):
    summary = open("3dmodels/summary.txt","a")
    
    summary.write("The parameter of the "+name+" orbit is a positive root of the following polynomial: "+str(polynomial)+"\n")
    
    summary.close()

def ExportAllFacetings(vertices, group, directory, name, ignoreTriangles = False, poly = None):
    facetings = FacetAll(vertices, group, ignoreTriangles)
    for i in range(len(facetings)):
        print(name+"."+str(i))
        ExportToOFF(vertices, facetings[i], group, directory, name+"."+str(i))
        
    if len(facetings) > 0:
        WriteSummary(name,poly)
        