#enumerate1D-A.py
#Part of the enumeration of noble polyhedra in 1D orbit types. 
#This file initializes the information needed for enumerate1D-B.wls
#and enumerate1D-C.py, as well as checking the typical equivalence
#classes for facetings. Together the 3 files complete the enumeration
#in this case.

print("Starting!")

from symbolic1D import *
import groups
import time
from faceting import FacetMinimalEquivalenceClass

startTime = time.time()

print("Obtaining volume configurations for small orbit types (tT,rT,rP,tO,tC,rC)...")
tTConf = VolumeConfiguration(tT)
rTConf = VolumeConfiguration(rT)
rPConf = VolumeConfiguration(rP)
tOConf = VolumeConfiguration(tO)
tCConf = VolumeConfiguration(tC)
rCConf = VolumeConfiguration(rC)

print("Obtaining volume configuration for tI orbit type...")
tIConf = VolumeConfiguration(tI)

print("Obtaining volume configuration for tD orbit type...")
tDConf = VolumeConfiguration(tD)

print("Obtaining volume configuration for rD orbit type...")
rDConf = VolumeConfiguration(rD)

confTime = time.time()
print("Volume configuration time:",confTime-startTime,"seconds.")

print("Faceting typical equivalence classes...")

minFacetings = []
for group in groups.tTGroups: #tT orbit type
    minFacetings += FacetMinimalEquivalenceClass(tT, MergePlanes(tTConf[0]), group)
for group in groups.rTGroups: #rT orbit type
    minFacetings += FacetMinimalEquivalenceClass(rT, MergePlanes(rTConf[0]), group)
for group in groups.rPGroups: #rP orbit type
    minFacetings += FacetMinimalEquivalenceClass(rP, MergePlanes(rPConf[0]), group)
for group in groups.tOGroups: #tO orbit type
    minFacetings += FacetMinimalEquivalenceClass(tO, MergePlanes(tOConf[0]), group)
for group in groups.tCGroups: #tC orbit type
    minFacetings += FacetMinimalEquivalenceClass(tC, MergePlanes(tCConf[0]), group)
for group in groups.rCGroups: #rC orbit type
    minFacetings += FacetMinimalEquivalenceClass(rC, MergePlanes(rCConf[0]), group)
for group in groups.tIGroups: #tI orbit type
    minFacetings += FacetMinimalEquivalenceClass(tI, MergePlanes(tIConf[0]), group)
for group in groups.tDGroups: #tD orbit type
    minFacetings += FacetMinimalEquivalenceClass(tD, MergePlanes(tDConf[0]), group)
for group in groups.rDGroups: #rD orbit type
    minFacetings += FacetMinimalEquivalenceClass(rD, MergePlanes(rDConf[0]), group)

#It turns out through these calculations that no noble polyhedra of
#this form exist, so we do not need to worry about exporting them.
if len(minFacetings) == 0:
    print("Found 0 facetings in minimum equivalence classes, as expected.")
else:
    raise Exception("Found unexpected faceting!")

minTime = time.time()
print("Typical faceting time:",minTime-confTime,"seconds.")

print("Exporting volume configurations...")
ExportConf(tTConf, "tT")
ExportConf(rTConf, "rT")
ExportConf(rPConf, "rP")
ExportConf(tOConf, "tO")
ExportConf(tCConf, "tC")
ExportConf(rCConf, "rC")
ExportConf(tIConf, "tI")
ExportConf(tDConf, "tD")
ExportConf(rDConf, "rD")
exportTime = time.time()

print("Configuration export time:",exportTime-minTime,"seconds.")
print("Finished. Total time taken:",exportTime-startTime,"seconds.")