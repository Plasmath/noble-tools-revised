from symbolic import *
import groups
import time

startTime = time.time()

#Volume configurations.
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
print("Total volume configuration time: %s seconds." % (confTime - startTime))

#Obtaining roots of entries within said volume configurations.
print("Obtaining parameter roots for small orbit types (tT,rT,rP,tO,tC,rC)...")
tTroots = Get1DOrbitRoots(tTConf, a)
rTroots = Get1DOrbitRoots(rTConf, a)
rProots = Get1DOrbitRoots(rPConf, a)
tOroots = Get1DOrbitRoots(tOConf, a)
tCroots = Get1DOrbitRoots(tCConf, a)
rCroots = Get1DOrbitRoots(rCConf, a)

print("Obtaining parameter roots for tI orbit type...")
tIroots = Get1DOrbitRoots(tIConf, a)

print("Obtaining parameter roots for tD orbit type...")
tDroots = Get1DOrbitRoots(tDConf, a)

print("Obtaining parameter roots for rD orbit type...")
rDroots = Get1DOrbitRoots(rDConf, a)

rootTime = time.time()
print("Total rootfinding time: %s seconds." % (rootTime - confTime))

#Enumerating facetings for all orbit types.
print("Enumerating small orbit types (tT,rT,rP,tO,tC,rC)...")
#tT orbit type
Enumerate1DOrbit(tTroots, tT, groups.tTGroup332, a, "3dmodels", "tT332")
Enumerate1DOrbit(tTroots, tT, groups.tTGroupStar332, a, "3dmodels", "tTStar332")
#rT orbit type
Enumerate1DOrbit(rTroots, rT, groups.rTGroup332, a, "3dmodels", "rT332")
Enumerate1DOrbit(rTroots, rT, groups.rTGroupStar332, a, "3dmodels", "rTStar332")
#rP orbit type
Enumerate1DOrbit(rProots, rP, groups.rPGroup3Star2, a, "3dmodels", "rP3Star2")
#tO orbit type
Enumerate1DOrbit(tOroots, tO, groups.tOGroup432, a, "3dmodels", "tO432")
Enumerate1DOrbit(tOroots, tO, groups.tOGroupStar432, a, "3dmodels", "tOStar432")
#tC orbit type
Enumerate1DOrbit(tCroots, tC, groups.tCGroup432, a, "3dmodels", "tC432")
Enumerate1DOrbit(tCroots, tC, groups.tCGroup3Star2, a, "3dmodels", "tC3Star2")
Enumerate1DOrbit(tCroots, tC, groups.tCGroupStar432, a, "3dmodels", "tCStar432")
#rC orbit type
Enumerate1DOrbit(rCroots, rC, groups.rCGroup432, a, "3dmodels", "rC432")
Enumerate1DOrbit(rCroots, rC, groups.rCGroup3Star2, a, "3dmodels", "rC3Star2")
Enumerate1DOrbit(rCroots, rC, groups.rCGroupStar432, a, "3dmodels", "rCStar432")

print("Enumerating tI orbit type...")
#tI orbit type
Enumerate1DOrbit(tIroots, tI, groups.tIGroup532, a, "3dmodels", "tI532")
Enumerate1DOrbit(tIroots, tI, groups.tIGroupStar532, a, "3dmodels", "tIStar532")

print("Enumerating tD orbit type...")
#tD orbit type
Enumerate1DOrbit(tDroots, tD, groups.tDGroup532, a, "3dmodels", "tD532")
Enumerate1DOrbit(tDroots, tD, groups.tDGroupStar532, a, "3dmodels", "tDStar532")

print("Enumerating rD orbit type...")
#tD orbit type
Enumerate1DOrbit(rDroots, rD, groups.rDGroup532, a, "3dmodels", "rD532")
Enumerate1DOrbit(rDroots, rD, groups.rDGroupStar532, a, "3dmodels", "rDStar532")

enumTime = time.time()
print("Total enumeration time: %s seconds." % (enumTime - rootTime))
print("Total execution time: %s seconds." % (enumTime - startTime))