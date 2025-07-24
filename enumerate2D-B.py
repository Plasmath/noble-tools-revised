#enumerate2D-B.py
#Find all critical pairs within coprime polynomials.

print("Starting!")

from symbolic2D import ImportCopr, MergeAllWithShared, FindCriticalPairs, ExportCriticalPairs
import time

startTime = time.time()

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