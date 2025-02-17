from faceting import Generate, FacetAll
from plot import MakePlot

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
    
    MakePlot(vertices, face, group, name)

def ExportAllFacetings(vertices, group, directory, name):
    facetings = FacetAll(vertices, group)
    for i in range(len(facetings)):
        ExportToOFF(vertices, facetings[i], group, directory, name+"."+str(i))
