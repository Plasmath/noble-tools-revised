print("Starting!")

from symbolic2D import *
import groups
import time

startTime = time.time()

#Volume configurations.
print("Obtaining volume configurations for small orbit types (sT,gT,gP,sC,gC)...")
sTConf = VolumeConfiguration(sT)
gTConf = VolumeConfiguration(gT)
gPConf = VolumeConfiguration(gP)
sCConf = VolumeConfiguration(sC)
gCConf = VolumeConfiguration(gC)

print("Obtaining volume configuration for sD orbit type...")
sDConf = VolumeConfiguration(sD)

print("Obtaining volume configuratrion for gD orbit type...")
gDConf = VolumeConfiguration(gD)

confTime = time.time()
print("Total volume configuration time: %s seconds." % (confTime - startTime))

#Coprime polynomials.
print("Obtaining coprime polynomials for small orbit types (sT,gT,gP,sC,gC)")
sTCopr = GetCopr(sT, sTConf)
gTCopr = GetCopr(gT, gTConf)
gPCopr = GetCopr(gP, gPConf)
sCCopr = GetCopr(sC, sCConf)
gCCopr = GetCopr(gC, gCConf)

print("Obtaining coprime polynomials for sD orbit type...")
sDCopr = GetCopr(sD, sDConf, extension = [sp.sqrt(2),sp.sqrt(5)])

print("Obtaining coprime polynomials for gD orbit type...")
gDCopr = GetCopr(gD, gDConf, extension = [sp.sqrt(2),sp.sqrt(5)])

coprTime = time.time()
print("Total coprime polynomial time: %s seconds." % (coprTime - confTime))

#Exporting
print("Exporting coprime polynomials...")
ExportCopr(sTCopr, "sT")
ExportCopr(gTCopr, "gT")
ExportCopr(gPCopr, "gP")
ExportCopr(sCCopr, "sC")
ExportCopr(gCCopr, "gC")
ExportCopr(sDCopr, "sD")
ExportCopr(gDCopr, "gD")
exportTime = time.time()

print("Total execution time: %s seconds." % (exportTime - startTime))