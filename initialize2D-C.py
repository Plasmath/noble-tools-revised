import sympy as sp
from sympy import sqrt
from symbolic2D import ImportCoprData, MergePlanes, MergeAllPlanes, ExportCopr
import config
import time

a = config.a
b = config.b

startTime = time.time()

#Reprocessing stored volume configurations
print("Importing volume configurations...")
confFile = open("data/conf.txt","r").read().split("\n")

sTConf = eval(confFile[0][5:])
gTConf = eval(confFile[1][5:])
gPConf = eval(confFile[2][5:])
sCConf = eval(confFile[3][5:])
gCConf = eval(confFile[4][5:])
sDConf = eval(confFile[5][5:])
gDConf = eval(confFile[6][5:])

importTimeA = time.time()
print("Volume configuration import time:",importTimeA - startTime)

#Importing ConfFactors and ConfFactorIndices files
print("Importing coprime polynomials...")

sTFactorData = ImportCoprData("sT", sTConf)
gTFactorData = ImportCoprData("gT", gTConf)
gPFactorData = ImportCoprData("gP", gPConf)
sCFactorData = ImportCoprData("sC", sCConf)
gCFactorData = ImportCoprData("gC", gCConf)
sDFactorData = ImportCoprData("sD", sDConf)
gDFactorData = ImportCoprData("gD", gDConf)

importTimeB = time.time()
print("Volume configuration import time:",importTimeB - importTimeA)

#Final plane merging and export
print("Finalizing coprime polynomials...")

sTCopr = MergeAllPlanes(sTFactorData)
gTCopr = MergeAllPlanes(gTFactorData)
gPCopr = MergeAllPlanes(gPFactorData)
sCCopr = MergeAllPlanes(sCFactorData)
gCCopr = MergeAllPlanes(gCFactorData)
sDCopr = MergeAllPlanes(sDFactorData)
gDCopr = MergeAllPlanes(gDFactorData)

finalizationTime = time.time()
print("Coprime finalization time:", finalizationTime - importTimeB)

print("Exporting...")
ExportCopr("sT",sTCopr, []) #sT has no planes with more than 3 vertices in the minimum equivalence class, so there are no shared planes.
ExportCopr("gT",gTCopr, MergePlanes(gTConf[0]))
ExportCopr("gP",gPCopr, MergePlanes(gPConf[0]))
ExportCopr("sC",sCCopr, MergePlanes(sCConf[0]))
ExportCopr("gC",gCCopr, MergePlanes(gCConf[0]))
ExportCopr("sD", sDCopr, MergePlanes(sDConf[0]))
ExportCopr("gD", gDCopr, MergePlanes(gDConf[0]))

endTime = time.time()
print("Total time:", endTime - startTime)