print("Starting!")
print("Note: Estimated runtime required to complete enumeration is approximately 30 minutes.")

from symbolic1D import *
import groups
import time

startTime = time.time()

#Volume configurations.
print("Obtaining coprime polynomials for small orbit types (tT,rT,rP,tO,tC,rC)...")
tTCopr = GetCopr(tT)
rTCopr = GetCopr(rT)
rPCopr = GetCopr(rP)
tOCopr = GetCopr(tO)
tCCopr = GetCopr(tC)
rCCopr = GetCopr(rC)

print("Obtaining coprime polynomials for tI orbit type...")
tICopr = GetCopr(tI, extension=[sp.sqrt(2),sp.sqrt(5)])

print("Obtaining coprime polynomials for tD orbit type...")
tDCopr = GetCopr(tD, extension=[sp.sqrt(2),sp.sqrt(5)])

print("Obtaining coprime polynomials for rD orbit type...")
rDCopr = GetCopr(rD, extension=[sp.sqrt(2),sp.sqrt(5)])

coprTime = time.time()
print("Total coprime polynomial time: %s seconds." % (coprTime - startTime))

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
print("Total candidate time: %s seconds." % (candTime - coprTime))

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
Export1DOrbitTypeFacetings(tI, tICandidatesStar432, groups.tIGroupStar532, "3dmodels", "tIStar532")
Export1DOrbitTypeFacetings(tI, tICandidates432, groups.tIGroup532, "3dmodels", "tI532")

print("Obtaining realizations for tD orbit type...")
Export1DOrbitTypeFacetings(tD, tDCandidatesStar432, groups.tDGroupStar532, "3dmodels", "tDStar532")
Export1DOrbitTypeFacetings(tD, tDCandidates432, groups.tDGroup532, "3dmodels", "tD532")

print("Obtaining realizations for rD orbit type...")
Export1DOrbitTypeFacetings(rD, rDCandidatesStar432, groups.rDGroupStar532, "3dmodels", "rDStar532")
Export1DOrbitTypeFacetings(rD, rDCandidates432, groups.rDGroup532, "3dmodels", "rD532")

enumTime = time.time()
print("Total realization time: %s seconds." % (enumTime - candTime))
print("Total execution time: %s seconds." % (enumTime - startTime))