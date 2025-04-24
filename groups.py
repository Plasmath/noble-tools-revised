#groups.py
#Code used to generate the given permutation groups for each orbit type.
#Note that the standard cycle notation is not used as other forms of storing the permutation are more computationally efficient.
#Reference Table 1 for more details on which orbit types are acted transitively on by which point groups.

from vertices import *
import numpy as np

def GetIndexOfImage(base,p,mat): #Get the index of mat*p in the list of points in the given base orbit.
    image = RoundArray(mat.dot(p))
    for i in range(len(base)):
        if np.all(RoundArray(base[i]) == image):
            return i

def GetPermutation(base,mat): #Get generating permutation given the generator matrix and base orbit.
    return tuple(GetIndexOfImage(base,p,mat) for p in base)

def Compose(perm1,perm2): #Compose two permutations to get a new permutation.
    return tuple(perm1[perm2[i]] for i in range(len(perm2)))

def GeneratePermutationGroup(base, generators, N): #Generate the permutation group for a given orbit type.
    
    perm1 = GetPermutation(base, generators[0])
    perm2 = GetPermutation(base, generators[1])
    perm3 = GetPermutation(base, generators[2])
    
    G = {perm1,perm2,perm3}
    
    for i in range(N):
        Gnew = G.copy()
        for g in G:
            Gnew.add(Compose(perm1,g))
            Gnew.add(Compose(perm2,g))
            Gnew.add(Compose(perm3,g))
        G = Gnew.copy()
    return G

#Generating matrices for the different symmetry groups.
genStar332 = [r1tet, r2tet, r3tet]
gen332 = [r1tet.dot(r2tet), r1tet.dot(r3tet), r2tet.dot(r3tet)]
genStar432 = [r1oct, r2oct, r3oct]
gen432 = [r1oct.dot(r2oct), r1oct.dot(r3oct), r2oct.dot(r3oct)]
genStar532 = [r1ico, r2ico, r3ico]
gen532 = [r1ico.dot(r2ico), r1ico.dot(r3ico), r2ico.dot(r3ico)]
gen3Star2 = [r1oct, r2oct.dot(r1oct.dot(r2oct)), r3oct.dot(r2oct)]

"""

0 Degrees of Freedom

"""

#T orbit type.
Tbase = GenerateVStar332(1, 0, 0)

TGroupStar332 = GeneratePermutationGroup(Tbase, genStar332, 10 )
TGroup332 = GeneratePermutationGroup(Tbase, gen332, 10)

#O orbit type.
Obase = GenerateVStar432(0, 0, 1)

OGroupStar332 = GeneratePermutationGroup(Obase, genStar332, 10 )
OGroup332 = GeneratePermutationGroup(Obase, gen332, 10 )
OGroup3Star2 = GeneratePermutationGroup(Obase, gen3Star2, 10 )
OGroupStar432 = GeneratePermutationGroup(Obase, genStar432, 10)
OGroup432 = GeneratePermutationGroup(Obase, gen432, 10)

#CO orbit type. The base orbit under the *332 subgroups is included within rT,
#so facetings under those symmetries are enumerated when rT is.
CObase = GenerateVStar432(0, 1, 0)

COGroup3Star2 = GeneratePermutationGroup(CObase, gen3Star2, 10)
COGroupStar432 = GeneratePermutationGroup(CObase, genStar432, 10)
COGroup432 = GeneratePermutationGroup(CObase, gen432, 10)

#C orbit type.
Cbase = GenerateVStar432(1, 0, 0)

CGroup3Star2 = GeneratePermutationGroup(Cbase, gen3Star2, 10)
CGroupStar432 = GeneratePermutationGroup(Cbase, genStar432, 10)
CGroup432 = GeneratePermutationGroup(Cbase, gen432, 10)

#I orbit type. The base orbit under the subgroups of *332 and 3*2 is included within rP,
#so facetings under those symmetries are enumerated when rP is.
Ibase = GenerateVStar532(0, 1, 0)

IGroupStar532 = GeneratePermutationGroup(Ibase, genStar532, 14)
IGroup532 = GeneratePermutationGroup(Ibase, gen532, 14)

#ID orbit type.
IDbase = GenerateVStar532(0, 0, 1)

IDGroupStar532 = GeneratePermutationGroup(IDbase, genStar532, 14)
IDGroup532 = GeneratePermutationGroup(IDbase, gen532, 14)

#D orbit type.
Dbase = GenerateVStar532(1, 0, 0)

DGroupStar532 = GeneratePermutationGroup(Dbase, genStar532, 18)
DGroup532 = GeneratePermutationGroup(Dbase, gen532, 14)

"""

1 Degree of Freedom

"""

#tT orbit type.
tTbase = GenerateVStar332(0, 1, 1)

tTGroupStar332 = GeneratePermutationGroup(tTbase, genStar332, 10)
tTGroup332 = GeneratePermutationGroup(tTbase, gen332, 10)
tTGroups = [tTGroupStar332,tTGroup332]

#rT orbit type.
rTbase = GenerateVStar332(1, 1, 0)

rTGroupStar332 = GeneratePermutationGroup(rTbase, genStar332, 10)
rTGroup332 = GeneratePermutationGroup(rTbase, gen332, 10)
rTGroups = [rTGroupStar332,rTGroup332]

#rP orbit type. The base orbit under 332 symmetry is included within sT,
#so facetings under that symmetry are enumerated when sT is.
rPbase = GenerateV3Star2(0, 1, 1)

rPGroup3Star2 = GeneratePermutationGroup(rPbase, gen3Star2, 10)
rPGroups = [rPGroup3Star2]

#tO orbit type. The base orbit under *332 symmetry is included within gT,
#so facetings under that symmetry are enumerated when gT is.
tObase = GenerateVStar432(0, 1, 1)

tOGroupStar432 = GeneratePermutationGroup(tObase, genStar432, 10)
tOGroup432 = GeneratePermutationGroup(tObase, gen432, 10)
tOGroups = [tOGroupStar432,tOGroup432]

#tC orbit type.
tCbase = GenerateVStar432(1, 1, 0)

tCGroup3Star2 = GeneratePermutationGroup(tCbase, gen3Star2, 10)
tCGroupStar432 = GeneratePermutationGroup(tCbase, genStar432, 10)
tCGroup432 = GeneratePermutationGroup(tCbase, gen432, 10)
tCGroups = [tCGroupStar432,tCGroup432,tCGroup3Star2]

#rC orbit type.
rCbase = GenerateVStar432(1, 0, 1)

rCGroup3Star2 = GeneratePermutationGroup(rCbase, gen3Star2, 10)
rCGroupStar432 = GeneratePermutationGroup(rCbase, genStar432, 10)
rCGroup432 = GeneratePermutationGroup(rCbase, gen432, 10)
rCGroups = [rCGroupStar432,rCGroup432,rCGroup3Star2]

#tI orbit type.
tIbase = GenerateVStar532(0, 1, 1)

tIGroupStar532 = GeneratePermutationGroup(tIbase, genStar532, 14)
tIGroup532 = GeneratePermutationGroup(tIbase, gen532, 14)
tIGroups = [tIGroupStar532,tIGroup532]

#tD orbit type.
tDbase = GenerateVStar532(1, 0, 1)

tDGroupStar532 = GeneratePermutationGroup(tDbase, genStar532, 14)
tDGroup532 = GeneratePermutationGroup(tDbase, gen532, 14)
tDGroups = [tDGroupStar532,tDGroup532]

#rD orbit type.
rDbase = GenerateVStar532(1, 1, 0)

rDGroupStar532 = GeneratePermutationGroup(rDbase, genStar532, 14)
rDGroup532 = GeneratePermutationGroup(rDbase, gen532, 14)
rDGroups = [rDGroupStar532,rDGroup532]

"""

2 Degrees of Freedom

"""

#sT orbit type.
sTbase = GenerateV332(1, 1, 1)

sTGroup332 = GeneratePermutationGroup(sTbase, gen332, 10)

#gT orbit type.
gTbase = GenerateVStar332(1, 1, 1)

gTGroupStar332 = GeneratePermutationGroup(gTbase, genStar332, 10)

#gP orbit type.
gPbase = GenerateV3Star2(1, 1, 1)

gPGroup3Star2 = GeneratePermutationGroup(gPbase, gen3Star2, 10)

#sC orbit type.
sCbase = GenerateV432(1, 1, 1)

sCGroup432 = GeneratePermutationGroup(sCbase, gen432, 10)

#gC orbit type.
gCbase = GenerateVStar432(1, 1, 1)

gCGroupStar432 = GeneratePermutationGroup(gCbase, genStar432, 10)

#sD orbit type.
sDbase = GenerateV532(1, 1, 1)

sDGroup532 = GeneratePermutationGroup(sDbase, gen532, 14)

#gD orbit type.
gDbase = GenerateVStar532(1, 1, 1)

gDGroupStar532 = GeneratePermutationGroup(gDbase, genStar532, 14)