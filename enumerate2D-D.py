#enumerate2D-D.py
#Converts the output of enumerate2D-C.wls to a more convenient format.

print("Starting!")

from symbolic2D import ImportCopr, MergeAllWithShared, ImportPairData, FindIsolatedPairs
import time

startTime = time.time()

print("Importing...")

sTPairData = ImportPairData("sT")
gTPairData = ImportPairData("gT")
gPPairData = ImportPairData("gP")
sCPairData = ImportPairData("sC")
gCPairData = ImportPairData("gC")
sDPairData = ImportPairData("sD")
gDPairData = ImportPairData("gD")

importTime = time.time()
print("Import time: %s seconds." % (importTime - startTime))

print("Importing pair data...")
