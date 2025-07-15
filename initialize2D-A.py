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

#Exporting
print("Exporting volume configurations...")
ExportConf(sTConf, "sT")
ExportConf(gTConf, "gT")
ExportConf(gPConf, "gP")
ExportConf(sCConf, "sC")
ExportConf(gCConf, "gC")
ExportConf(sDConf, "sD")
ExportConf(gDConf, "gD")
exportTime = time.time()

print("Total execution time: %s seconds." % (exportTime - startTime))