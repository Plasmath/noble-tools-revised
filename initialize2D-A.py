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
ExportConf(sTConf, "sT", dim = 2)
ExportConf(gTConf, "gT", dim = 2)
ExportConf(gPConf, "gP", dim = 2)
ExportConf(sCConf, "sC", dim = 2)
ExportConf(gCConf, "gC", dim = 2)
ExportConf(sDConf, "sD", dim = 2)
ExportConf(gDConf, "gD", dim = 2)
exportTime = time.time()

print("Total execution time: %s seconds." % (exportTime - startTime))