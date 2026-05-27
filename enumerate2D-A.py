#enumerate2D-A.py
#Enumerates noble polyhedra in orbit types with 2 degrees of freedom,
#but not within a maximal equivalence class.

print("Starting!")

from sympy import sqrt
from config import a,b
import groups
import time
import symbolic
from symbolic2D import Get2DCoprCandidates, ImportCopr
from faceting import FacetMinimalEquivalenceClass

print("Importing...")

startTime = time.time()

sTCopr, sTShared = ImportCopr("sT")
gTCopr, gTShared = ImportCopr("gT")
gPCopr, gPShared = ImportCopr("gP")
sCCopr, sCShared = ImportCopr("sC")
gCCopr, gCShared = ImportCopr("gC")
sDCopr, sDShared = ImportCopr("sD")
gDCopr, gDShared = ImportCopr("gD")

importTime = time.time()
print("Import time: %s seconds." % (importTime - startTime))

#Check the typical equivalence classes for facetings
print("Checking typical equivalence classes...")
minFacetings = []
minFacetings += FacetMinimalEquivalenceClass(symbolic.sT, sTShared, groups.sTGroup332) #sT orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.gT, gTShared, groups.gTGroupStar332) #gT orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.gP, gPShared, groups.gPGroup3Star2) #gP orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.sC, sCShared, groups.sCGroup432) #sC orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.gC, gCShared, groups.gCGroupStar432) #gC orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.sD, sDShared, groups.sDGroup532) #sD orbit type
minFacetings += FacetMinimalEquivalenceClass(symbolic.gD, gDShared, groups.gDGroupStar532) #gD orbit type

#Like in the 1D case, no noble polyhedra of this form turn out to exist.
#Therefore, we do not need to worry about exporting them, but we have this
#error to handle the case where any are found.
if len(minFacetings) == 0:
    print("Found 0 facetings in minimum equivalence classes, as expected.")
else:
    raise Exception("Found unexpected faceting!")

minEquivTime = time.time()
print("Minimum equivalence class faceting time: %s seconds." % (minEquivTime - importTime))

#Check the other nonmaximal equivalence classes for facetings - there exists exactly
#one of these for each coprime polynomial of each orbit type.
print("Obtaining candidates for other nonmaximal equivalence classes...")

candidates = []
candidates += Get2DCoprCandidates(symbolic.sT, sTCopr, sTShared, groups.sTGroup332)
candidates += Get2DCoprCandidates(symbolic.gT, gTCopr, gTShared, groups.gTGroupStar332)
candidates += Get2DCoprCandidates(symbolic.gP, gPCopr, gPShared, groups.gPGroup3Star2)
candidates += Get2DCoprCandidates(symbolic.sC, sCCopr, sCShared, groups.sCGroup432)
candidates += Get2DCoprCandidates(symbolic.gC, gCCopr, gCShared, groups.gCGroupStar432)
candidates += Get2DCoprCandidates(symbolic.sD, sDCopr, sDShared, groups.sDGroup532)
candidates += Get2DCoprCandidates(symbolic.gD, gDCopr, gDShared, groups.gDGroupStar532)

#Again, no noble polyhedra of this form exist, but we have this exception here for logical soundness.
if len(candidates) == 0:
    print("Found 0 candidates in other nonmaximal equivalence classes, as expected.")
else:
    raise Exception("Found undexpected faceting!")

candidateTime = time.time()
print("Candidate time:",candidateTime - minEquivTime)

print("Execution successful: no noble polyhedra found in these equivalence classes.")
print("Total time: %s seconds." % (candidateTime - startTime))