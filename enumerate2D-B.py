#Find all critical pairs within coprime polynomials.

from symbolic2D import ImportCopr, FindCriticalPairs
import time

startTime = time.time()

print("Starting!")

sTCopr, sTShared = ImportCopr("sT")
gTCopr, gTShared = ImportCopr("gT")
gPCopr, gPShared = ImportCopr("gP")
sCCopr, sCShared = ImportCopr("sC")
gCCopr, gCShared = ImportCopr("gC")
sDCopr, sDShared = ImportCopr("sD")
gDCopr, gDShared = ImportCopr("gD")

importTime = time.time()
print("Import time:",importTime - startTime)

x = lambda pairs, copr : print(len(pairs), len(pairs)*2/(len(copr)*(len(copr)-1)))

sTPairs = FindCriticalPairs(sTCopr, sTShared)
x(sTPairs,sTCopr)

gTPairs = FindCriticalPairs(gTCopr, gTShared)
x(gTPairs,gTCopr)

gPPairs = FindCriticalPairs(gPCopr, gPShared)
x(gPPairs,gPCopr)

sCPairs = FindCriticalPairs(sCCopr, sCShared)
x(sCPairs,sCCopr)

gCPairs = FindCriticalPairs(gCCopr, gCShared)
x(gCPairs,gCCopr)

sDPairs = FindCriticalPairs(sDCopr, sDShared)
x(sDPairs,sDCopr)

gDPairs = FindCriticalPairs(gDCopr, gDShared)
x(gDPairs,gDCopr)
