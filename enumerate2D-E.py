# enumerate2D-E.py
# The final step where we take the abstract faceting data found by previous
# files and check to see which of these leads to real polyhedra. This data
# is then exported to the 3dmodels folder.

from symbolic2D import ImportAtlasData, ImportIntersectionData, ImportFacetingData, DetermineIntersectionData, Get2DOrbitTypeFacetings, Export2DOrbitTypeFacetings
import symbolic2D
import time
import groups

print("Starting!")
startTime = time.time()

# Obtain the atlases of roots for the orbit types
sTAtlasData = ImportAtlasData("sT")
gTAtlasData = ImportAtlasData("gT")
gPAtlasData = ImportAtlasData("gP")
sCAtlasData = ImportAtlasData("sC")
gCAtlasData = ImportAtlasData("gC")
sDAtlasData = ImportAtlasData("sD")
gDAtlasData = ImportAtlasData("gD")

sTIntersections = ImportIntersectionData("sT")
gTIntersections = ImportIntersectionData("gT")
gPIntersections = ImportIntersectionData("gP")
sCIntersections = ImportIntersectionData("sC")
gCIntersections = ImportIntersectionData("gC")
sDIntersections = ImportIntersectionData("sD")
gDIntersections = ImportIntersectionData("gD")

sTPoly, sTCopr, sTSharedPlanes = ImportFacetingData("sT")
gTPoly, gTCopr, gTSharedPlanes = ImportFacetingData("gT")
gPPoly, gPCopr, gPSharedPlanes = ImportFacetingData("gP")
sCPoly, sCCopr, sCSharedPlanes = ImportFacetingData("sC")
gCPoly, gCCopr, gCSharedPlanes = ImportFacetingData("gC")
sDPoly, sDCopr, sDSharedPlanes = ImportFacetingData("sD")
gDPoly, gDCopr, gDSharedPlanes = ImportFacetingData("gD")

importTime = time.time()
print("Import time:",importTime - startTime)

#Determine planes for faceting
sTIntersectionData = DetermineIntersectionData(sTIntersections, sTPoly, sTCopr, sTSharedPlanes)
gTIntersectionData = DetermineIntersectionData(gTIntersections, gTPoly, gTCopr, gTSharedPlanes)
gPIntersectionData = DetermineIntersectionData(gPIntersections, gPPoly, gPCopr, gPSharedPlanes)
sCIntersectionData = DetermineIntersectionData(sCIntersections, sCPoly, sCCopr, sCSharedPlanes)
gCIntersectionData = DetermineIntersectionData(gCIntersections, gCPoly, gCCopr, gCSharedPlanes)
sDIntersectionData = DetermineIntersectionData(sDIntersections, sDPoly, sDCopr, sDSharedPlanes)
gDIntersectionData = DetermineIntersectionData(gDIntersections, gDPoly, gDCopr, gDSharedPlanes)

planeTime = time.time()
print("Time spent preparing for faceting:", planeTime - importTime)

#Determine facetings
sTFacetings = Get2DOrbitTypeFacetings(symbolic2D.sT, sTIntersectionData, groups.sTGroup332)
gTFacetings = Get2DOrbitTypeFacetings(symbolic2D.gT, gTIntersectionData, groups.gTGroupStar332)
gPFacetings = Get2DOrbitTypeFacetings(symbolic2D.gP, gPIntersectionData, groups.gPGroup3Star2)
sCFacetings = Get2DOrbitTypeFacetings(symbolic2D.sC, sCIntersectionData, groups.sCGroup432)
gCFacetings = Get2DOrbitTypeFacetings(symbolic2D.gC, gCIntersectionData, groups.gCGroupStar432)
sDFacetings = Get2DOrbitTypeFacetings(symbolic2D.sD, sDIntersectionData, groups.sDGroup532)
gDFacetings = Get2DOrbitTypeFacetings(symbolic2D.gD, gDIntersectionData, groups.gDGroupStar532)

facetingTime = time.time()
print("Time spent faceting:", facetingTime - planeTime)

#Export facetings
Export2DOrbitTypeFacetings(symbolic2D.sT, sTAtlasData, sTFacetings, groups.sTGroup332, "3dmodels", "sT")
Export2DOrbitTypeFacetings(symbolic2D.gT, gTAtlasData, gTFacetings, groups.gTGroupStar332, "3dmodels", "gT")
Export2DOrbitTypeFacetings(symbolic2D.gP, gPAtlasData, gPFacetings, groups.gPGroup3Star2, "3dmodels", "gP")
Export2DOrbitTypeFacetings(symbolic2D.sC, sCAtlasData, sCFacetings, groups.sCGroup432, "3dmodels", "sC")
Export2DOrbitTypeFacetings(symbolic2D.gC, gCAtlasData, gCFacetings, groups.gCGroupStar432, "3dmodels", "gC")
Export2DOrbitTypeFacetings(symbolic2D.sD, sDAtlasData, sDFacetings, groups.sDGroup532, "3dmodels", "sD")
Export2DOrbitTypeFacetings(symbolic2D.gD, gDAtlasData, gDFacetings, groups.gDGroupStar532, "3dmodels", "gD")