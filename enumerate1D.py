from symbolic import *
import groups

print("Enumerating tT orbit type...")
#tT orbit type
Enumerate1DOrbit(tT, groups.tTGroup332, a, "3dmodels", "tT332")
Enumerate1DOrbit(tT, groups.tTGroupStar332, a, "3dmodels", "tTStar332")

print("Enumerating rT orbit type...")
#rT orbit type
Enumerate1DOrbit(rT, groups.rTGroup332, a, "3dmodels", "rT332")
Enumerate1DOrbit(rT, groups.rTGroupStar332, a, "3dmodels", "rTStar332")

print("Enumerating rP orbit type...")
#rP orbit type
Enumerate1DOrbit(rP, groups.rPGroup3Star2, a, "3dmodels", "rP3Star2")

print("Enumerating tO orbit type...")
#tO orbit type
Enumerate1DOrbit(tO, groups.tOGroup432, a, "3dmodels", "tO432")
Enumerate1DOrbit(tO, groups.tOGroupStar432, a, "3dmodels", "tOStar432")

print("Enumerating tC orbit type...")
#tC orbit type
Enumerate1DOrbit(tC, groups.tCGroup432, a, "3dmodels", "tC432")
Enumerate1DOrbit(tC, groups.tCGroup3Star2, a, "3dmodels", "tC3Star2")
Enumerate1DOrbit(tC, groups.tCGroupStar432, a, "3dmodels", "tCStar432")

print("Enumerating rC orbit type...")
#rC orbit type
Enumerate1DOrbit(rC, groups.rCGroup432, a, "3dmodels", "rC432")
Enumerate1DOrbit(rC, groups.rCGroup3Star2, a, "3dmodels", "rC3Star2")
Enumerate1DOrbit(rC, groups.rCGroupStar432, a, "3dmodels", "rCStar432")

print("Enumerating tI orbit type...")
#tI orbit type
Enumerate1DOrbit(tI, groups.tIGroup532, a, "3dmodels", "tI532")
Enumerate1DOrbit(tI, groups.tIGroupStar532, a, "3dmodels", "tIStar532")

print("Enumerating tD orbit type")
#tD orbit type
Enumerate1DOrbit(tD, groups.tDGroup532, a, "3dmodels", "tD532")
Enumerate1DOrbit(tD, groups.tDGroupStar532, a, "3dmodels", "tDStar532")

print("Enumerating rD orbit type")
#tD orbit type
Enumerate1DOrbit(rD, groups.rDGroup532, a, "3dmodels", "rD532")
Enumerate1DOrbit(rD, groups.rDGroupStar532, a, "3dmodels", "rDStar532")