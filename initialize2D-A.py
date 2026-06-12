#initialize2D-A.py
#The first of a group of three files which obtain the coprime polynomials
#for the orbit types with 2 degrees of freedom.

print("Starting!")

from symbolic2D import *
import groups
import time

startTime = time.time()

#Generate the different volume configurations.
print("Obtaining volume configurations for small orbit types (sT,gT,gP,sC,gC)...")
sTConf = VolumeConfiguration(sT) #sT orbit type
gTConf = VolumeConfiguration(gT) #gT orbit type
gPConf = VolumeConfiguration(gP) #gP orbit type
sCConf = VolumeConfiguration(sC) #sC orbit type
gCConf = VolumeConfiguration(gC) #gC orbit type

print("Obtaining volume configuration for sD orbit type...")
sDConf = VolumeConfiguration(sD) #sD orbit type

print("Obtaining volume configuratrion for gD orbit type...")
gDConf = VolumeConfiguration(gD) #gD orbit tyoe

confTime = time.time()
print("Total volume configuration time: %s seconds." % (confTime - startTime))

#Export the volume configurations that we just calculated.
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