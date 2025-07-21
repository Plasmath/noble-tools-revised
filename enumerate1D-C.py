#enumerate1D-C.py
#Part of the enumeration of noble polyhedra in 1D orbit types. 
#This file completes the enumeration given the information from
#enumerate1D-B.wls and enumerate1D-A.py.

print("Starting!")

from symbolic1D import *
import groups
import time

startTime = time.time()

#Reprocessing stored volume configurations
print("Importing volume configurations...")
confFile = open("data/conf1d.txt","r").read().split("\n")

tTConf = eval(confFile[0][5:])
rTConf = eval(confFile[1][5:])
rPConf = eval(confFile[2][5:])
tOConf = eval(confFile[3][5:])
tCConf = eval(confFile[4][5:])
rCConf = eval(confFile[5][5:])
tIConf = eval(confFile[6][5:])
tDConf = eval(confFile[7][5:])
rDConf = eval(confFile[8][5:])

importTimeA = time.time()
print("Volume configuration import time: %s seconds." % (importTimeA - startTime))

#Importing ConfFactors and ConfFactorIndices files
print("Importing coprime polynomials...")

tTFactorData = ImportCoprData("tT", tTConf)
rTFactorData = ImportCoprData("rT", rTConf)
rPFactorData = ImportCoprData("rP", rPConf)
tOFactorData = ImportCoprData("tO", tOConf)
tCFactorData = ImportCoprData("tC", tCConf)
rCFactorData = ImportCoprData("rC", rCConf)
tIFactorData = ImportCoprData("tI", tIConf)
tDFactorData = ImportCoprData("tD", tDConf)
rDFactorData = ImportCoprData("rD", rDConf)

importTimeB = time.time()
print("Coprime polynomial import time: %s seconds." % (importTimeB - importTimeA))

#Final plane merging and export
print("Finalizing coprime polynomials...")

tTCopr = MergeAllPlanes(tTFactorData)
rTCopr = MergeAllPlanes(rTFactorData)
rPCopr = MergeAllPlanes(rPFactorData)
tOCopr = MergeAllPlanes(tOFactorData)
tCCopr = MergeAllPlanes(tCFactorData)
rCCopr = MergeAllPlanes(rCFactorData)
tICopr = MergeAllPlanes(tIFactorData)
tDCopr = MergeAllPlanes(tDFactorData)
rDCopr = MergeAllPlanes(rDFactorData)

finalizationTime = time.time()
print("Coprime finalization time: %s seconds." % (finalizationTime - importTimeB))

#Obtaining orbit candidates
print("Obtaining orbit candidates for small orbit types (tT,rT,rP,tO,tC,rC)...")
tTCandidatesStar332 = Get1DOrbitCandidates(tT, tTCopr, groups.tTGroupStar332)
tTCandidates332 = Get1DOrbitCandidates(tT, tTCopr, groups.tTGroup332)

rTCandidatesStar332 = Get1DOrbitCandidates(rT, rTCopr, groups.rTGroupStar332)
rTCandidates332 = Get1DOrbitCandidates(rT, rTCopr, groups.rTGroup332)

rPCandidates3Star2 = Get1DOrbitCandidates(rP, rPCopr, groups.rPGroup3Star2)

tOCandidatesStar432 = Get1DOrbitCandidates(tO, tOCopr, groups.tOGroupStar432)
tOCandidates432 = Get1DOrbitCandidates(tO, tOCopr, groups.tOGroup432)

tCCandidatesStar432 = Get1DOrbitCandidates(tC, tCCopr, groups.tCGroupStar432)
tCCandidates432 = Get1DOrbitCandidates(tC, tCCopr, groups.tCGroup432)
tCCandidates3Star2 = Get1DOrbitCandidates(tC, tCCopr, groups.tCGroup3Star2)

rCCandidatesStar432 = Get1DOrbitCandidates(rC, rCCopr, groups.rCGroupStar432)
rCCandidates432 = Get1DOrbitCandidates(rC, rCCopr, groups.rCGroup432)
rCCandidates3Star2 = Get1DOrbitCandidates(rC, rCCopr, groups.rCGroup3Star2)

print("Obtaining orbit candidates for tI orbit type...")
tICandidatesStar532 = Get1DOrbitCandidates(tI, tICopr, groups.tIGroupStar532)
tICandidates532 = Get1DOrbitCandidates(tI, tICopr, groups.tIGroup532)

print("Obtaining orbit candidates for tD orbit type...")
tDCandidatesStar532 = Get1DOrbitCandidates(tD, tDCopr, groups.tDGroupStar532)
tDCandidates532 = Get1DOrbitCandidates(tD, tDCopr, groups.tDGroup532)

print("Obtaining orbit candidates for rD orbit type...")
rDCandidatesStar532 = Get1DOrbitCandidates(rD, rDCopr, groups.rDGroupStar532)
rDCandidates532 = Get1DOrbitCandidates(rD, rDCopr, groups.rDGroup532)

candTime = time.time()
print("Total candidate time: %s seconds." % (candTime - finalizationTime))

#Obtaining realizations
print("Obtaining realizations for small orbit types (tT,rT,rP,tO,tC,rC)...")
Export1DOrbitTypeFacetings(tT, tTCandidatesStar332, groups.tTGroupStar332, "3dmodels", "tTStar332")
Export1DOrbitTypeFacetings(tT, tTCandidates332, groups.tTGroup332, "3dmodels", "tT332")

Export1DOrbitTypeFacetings(rT, rTCandidatesStar332, groups.rTGroupStar332, "3dmodels", "rTStar332")
Export1DOrbitTypeFacetings(rT, rTCandidates332, groups.rTGroup332, "3dmodels", "rT332")

Export1DOrbitTypeFacetings(rP, rPCandidates3Star2, groups.rPGroup3Star2, "3dmodels", "rP3Star2")

Export1DOrbitTypeFacetings(tO, tOCandidatesStar432, groups.tOGroupStar432, "3dmodels", "tOStar432")
Export1DOrbitTypeFacetings(tO, tOCandidates432, groups.tOGroup432, "3dmodels", "tO432")

Export1DOrbitTypeFacetings(tC, tCCandidatesStar432, groups.tCGroupStar432, "3dmodels", "tCStar432")
Export1DOrbitTypeFacetings(tC, tCCandidates432, groups.tCGroup432, "3dmodels", "tC432")
Export1DOrbitTypeFacetings(tC, tCCandidates3Star2, groups.tCGroup3Star2, "3dmodels", "tC3Star2")

Export1DOrbitTypeFacetings(rC, rCCandidatesStar432, groups.rCGroupStar432, "3dmodels", "rCStar432")
Export1DOrbitTypeFacetings(rC, rCCandidates432, groups.rCGroup432, "3dmodels", "rC432")
Export1DOrbitTypeFacetings(rC, rCCandidates3Star2, groups.rCGroup3Star2, "3dmodels", "rC3Star2")

print("Obtaining realizations for tI orbit type...")
Export1DOrbitTypeFacetings(tI, tICandidatesStar532, groups.tIGroupStar532, "3dmodels", "tIStar532")
Export1DOrbitTypeFacetings(tI, tICandidates532, groups.tIGroup532, "3dmodels", "tI532")

print("Obtaining realizations for tD orbit type...")
Export1DOrbitTypeFacetings(tD, tDCandidatesStar532, groups.tDGroupStar532, "3dmodels", "tDStar532")
Export1DOrbitTypeFacetings(tD, tDCandidates532, groups.tDGroup532, "3dmodels", "tD532")

print("Obtaining realizations for rD orbit type...")
Export1DOrbitTypeFacetings(rD, rDCandidatesStar532, groups.rDGroupStar532, "3dmodels", "rDStar532")
Export1DOrbitTypeFacetings(rD, rDCandidates532, groups.rDGroup532, "3dmodels", "rD532")

enumTime = time.time()
print("Total realization time: %s seconds." % (enumTime - candTime))
print("Total execution time: %s seconds." % (enumTime - startTime))