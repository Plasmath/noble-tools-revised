#Enumerate noble polyhedra in orbit types with 2 degrees of freedom,
#but not a maximal equivalence class.
from sympy import sqrt
from config import a,b
import groups
import time
import symbolic
from symbolic2D import Get2DCoprCandidates, ImportCopr
from faceting import FacetMinimalEquivalenceClass

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

#Check the minimum equivalence classes for facetings
print("Checking the minimum equivalence classes...")

minFacetings = []
minFacetings += FacetMinimalEquivalenceClass(symbolic.sT, sTShared, groups.sTGroup332)
minFacetings += FacetMinimalEquivalenceClass(symbolic.gT, gTShared, groups.gTGroupStar332)
minFacetings += FacetMinimalEquivalenceClass(symbolic.gP, gPShared, groups.gPGroup3Star2)
minFacetings += FacetMinimalEquivalenceClass(symbolic.sC, sCShared, groups.sCGroup432)
minFacetings += FacetMinimalEquivalenceClass(symbolic.gC, gCShared, groups.gCGroupStar432)
minFacetings += FacetMinimalEquivalenceClass(symbolic.sD, sDShared, groups.sDGroup532)
minFacetings += FacetMinimalEquivalenceClass(symbolic.gD, gDShared, groups.gDGroupStar532)

#It turns out through these calculations that no noble polyhedra of
#this form exist, so we do not need to worry about exporting them.
if len(minFacetings) == 0:
    print("Found 0 facetings in minimum equivalence classes, as expected.")
else:
    raise Exception("Found unexpected faceting!")

minEquivTime = time.time()
print("Minimum equivalence class faceting time:",minEquivTime - importTime)

#Check the other nonmaximal equivalence classes for facetings
print("Obtaining candidates for other nonmaximal equivalence classes...")

candidates = []
candidates += Get2DCoprCandidates(symbolic.sT, sTCopr, sTShared, groups.sTGroup332)
candidates += Get2DCoprCandidates(symbolic.gT, gTCopr, gTShared, groups.gTGroupStar332)
candidates += Get2DCoprCandidates(symbolic.gP, gPCopr, gPShared, groups.gPGroup3Star2)
candidates += Get2DCoprCandidates(symbolic.sC, sCCopr, sCShared, groups.sCGroup432)
candidates += Get2DCoprCandidates(symbolic.gC, gCCopr, gCShared, groups.gCGroupStar432)
candidates += Get2DCoprCandidates(symbolic.sD, sDCopr, sDShared, groups.sDGroup532)
candidates += Get2DCoprCandidates(symbolic.gD, gDCopr, gDShared, groups.gDGroupStar532)

#It turns out through these calculations that no noble polyhedra of this form
#exist, so we do not need to worry about determining if they can be realized.
if len(candidates) == 0:
    print("Found 0 candidates in other nonmaximal equivalence classes, as expected.")
else:
    raise Exception("Found undexpected faceting!")

candidateTime = time.time()
print("Candidate time:",candidateTime - minEquivTime)

print("Execution successful: no noble polyhedra found in these equivalence classes.")
print("Total time:", candidateTime - startTime)