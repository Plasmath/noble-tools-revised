#enumerate0D.py
#Facets all orbit types with 0 degrees of freedom and exports
#their facetings as OFF files. This will contain some duplicates
#as some facetings are noble under multiple symmetries.

print("Starting!")

import groups
import time
from export import ExportAllFacetings

startTime = time.time()
print("Faceting orbits...")

T = groups.Tbase
O = groups.Obase
C = groups.Cbase
CO = groups.CObase
D = groups.Dbase
I = groups.Ibase
ID = groups.IDbase

#T orbit type
ExportAllFacetings(T, groups.TGroupStar332, "3dmodels", "TStar332")
ExportAllFacetings(T, groups.TGroup332, "3dmodels", "T332")

#O orbit type
ExportAllFacetings(O, groups.OGroupStar332, "3dmodels", "OStar332")
ExportAllFacetings(O, groups.OGroup332, "3dmodels", "O332")
ExportAllFacetings(O, groups.OGroup3Star2, "3dmodels", "O3Star2")
ExportAllFacetings(O, groups.OGroupStar432, "3dmodels", "OStar432")
ExportAllFacetings(O, groups.OGroup432, "3dmodels", "O432")

#C orbit type
ExportAllFacetings(C, groups.CGroup3Star2, "3dmodels", "C3Star2")
ExportAllFacetings(C, groups.CGroupStar432, "3dmodels", "CStar432")
ExportAllFacetings(C, groups.CGroup432, "3dmodels", "C432")

#CO orbit type
ExportAllFacetings(CO, groups.COGroup3Star2, "3dmodels", "CO3Star2")
ExportAllFacetings(CO, groups.COGroupStar432, "3dmodels", "COStar432")

ExportAllFacetings(CO, groups.COGroup432, "3dmodels", "CO432")

#I orbit type
ExportAllFacetings(I, groups.IGroupStar532, "3dmodels", "IStar532")
ExportAllFacetings(I, groups.IGroup532, "3dmodels", "I532")

#D orbit type
ExportAllFacetings(D, groups.DGroupStar532, "3dmodels", "DStar532")
ExportAllFacetings(D, groups.DGroup532, "3dmodels", "D532")

#ID orbit type
ExportAllFacetings(ID, groups.IDGroupStar532, "3dmodels", "IDStar532")
ExportAllFacetings(ID, groups.IDGroup532, "3dmodels", "ID532")

print("Finished.")
print("Total time taken:",time.time()-startTime)