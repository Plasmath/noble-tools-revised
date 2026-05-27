#enumerate2D-B.py
#Find pairs of coprime polynomials which may lead to new planes when both critical curves overlap.

print("Starting!")

from symbolic2D import ImportCopr, MergeAllWithShared, FindCriticalPairs, ExportCriticalPairs
import time

startTime = time.time()

#Import coprime polynomials for the orbit types and obtain the planes for their faceting data.
print("Importing...")
sTCopr, sTShared = ImportCopr("sT")
gTCopr, gTShared = ImportCopr("gT")
gPCopr, gPShared = ImportCopr("gP")
sCCopr, sCShared = ImportCopr("sC")
gCCopr, gCShared = ImportCopr("gC")
sDCopr, sDShared = ImportCopr("sD")
gDCopr, gDShared = ImportCopr("gD")
sTCopr = MergeAllWithShared(sTCopr, sTShared)
gTCopr = MergeAllWithShared(gTCopr, gTShared)
gPCopr = MergeAllWithShared(gPCopr, gPShared)
sCCopr = MergeAllWithShared(sCCopr, sCShared)
gCCopr = MergeAllWithShared(gCCopr, gCShared)
sDCopr = MergeAllWithShared(sDCopr, sDShared)
gDCopr = MergeAllWithShared(gDCopr, gDShared)
importTime = time.time()
print("Import time: %s seconds." % (importTime - startTime))

#Find pairs of coprime polynomials such that larger planes will exist in orbits at the
#intersections of these two critical curves. We don't need to consider other intersections
#as they can't lead to new polygon options for faces and thus new noble polyhedra.
print("Obtaining critical pairs...")
sTPairs = FindCriticalPairs(sTCopr)
gTPairs = FindCriticalPairs(gTCopr)
gPPairs = FindCriticalPairs(gPCopr)
sCPairs = FindCriticalPairs(sCCopr)
gCPairs = FindCriticalPairs(gCCopr)
sDPairs = FindCriticalPairs(sDCopr)
gDPairs = FindCriticalPairs(gDCopr)
pairTime = time.time()
print("Pair finding time: %s seconds." % (pairTime - importTime))

#Export the critical pair data for each orbit type
print("Exporting pairs...")
ExportCriticalPairs("sT",sTPairs,list(sTCopr.keys()))
ExportCriticalPairs("gT",gTPairs,list(gTCopr.keys()))
ExportCriticalPairs("gP",gPPairs,list(gPCopr.keys()))
ExportCriticalPairs("sC",sCPairs,list(sCCopr.keys()))
ExportCriticalPairs("gC",gCPairs,list(gCCopr.keys()))
ExportCriticalPairs("sD",sDPairs,list(sDCopr.keys()))
ExportCriticalPairs("gD",gDPairs,list(gDCopr.keys()))

exportTime = time.time()
print("Export time: %s seconds." % (exportTime - pairTime))