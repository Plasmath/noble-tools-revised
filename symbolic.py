#Code that constructs the symbolic forms of the orbit types.
#This is only needed for the orbit types with one or more degrees of freedom
#in order to obtain exact minimal polynomials for their parameters.
import vertices
import config

config.useSymbolic = True

vertices.InitializeSymbolic()

#Parameters for the orbit types.
a = config.a
b = config.b

#Orbit types with 1 degree of freedom.
tT = vertices.GenerateVStar332(0, a, 1)
rT = vertices.GenerateVStar332(a, 1, 0)
rP = vertices.GenerateV3Star2 (0, a, 1)
tO = vertices.GenerateVStar432(0, a, 1)
tC = vertices.GenerateVStar432(a, 1, 0)
rC = vertices.GenerateVStar432(a, 0, 1)
tI = vertices.GenerateVStar532(0, a, 1)
tD = vertices.GenerateVStar532(a, 0, 1)
rD = vertices.GenerateVStar532(a, 1, 0)

#Orbit types with 2 degrees of freedom.
sT = vertices.GenerateV332(a, b, 1)
gT = vertices.GenerateVStar332(a, b, 1)
gP = vertices.GenerateV3Star2(a, b, 1)
sC = vertices.GenerateV432(a, b, 1)
gC = vertices.GenerateVStar432(a, b, 1)
sD = vertices.GenerateV532(a, b, 1)
gD = vertices.GenerateVStar532(a, b, 1)