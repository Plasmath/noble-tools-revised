#!/usr/bin/env wolframscript

Off[General::stop]

Print["Importing..."]

(* Import polynomial data from enumerate2D-B.py *)

polysT = Import["data/sTCoprList.txt", "Table"]
polygT = Import["data/gTCoprList.txt", "Table"]
polygP = Import["data/gPCoprList.txt", "Table"]
polysC = Import["data/sCCoprList.txt", "Table"]
polygC = Import["data/gCCoprList.txt", "Table"]
polysD = Import["data/sDCoprList.txt", "Table"]
polygD = Import["data/gDCoprList.txt", "Table"]

pairsT = Import["data/sTPairs.txt", "Table"]
pairgT = Import["data/gTPairs.txt", "Table"]
pairgP = Import["data/gPPairs.txt", "Table"]
pairsC = Import["data/sCPairs.txt", "Table"]
pairgC = Import["data/gCPairs.txt", "Table"]
pairsD = Import["data/sDPairs.txt", "Table"]
pairgD = Import["data/gDPairs.txt", "Table"]

(* Import atlas data from enumerate2D-C.wls *)

sTAtlasData = ReadList["data/sTAtlas.txt"]
gTAtlasData = ReadList["data/gTAtlas.txt"]
gPAtlasData = ReadList["data/gPAtlas.txt"]
sCAtlasData = ReadList["data/sCAtlas.txt"]
gCAtlasData = ReadList["data/gCAtlas.txt"]
sDAtlasData = ReadList["data/sDAtlas.txt"]
gDAtlasData = ReadList["data/gDAtlas.txt"]

Print["Processing pair data..."]

GetRoot[poly_, atlaxApprox_, var_] := (
    For[ n = 1, n < PolynomialDegree[poly, var]+1, n++,
        If[ atlaxApprox == Root[poly,n] ]
    ];
)

FindIntersections[pairs_, poly_, atlasData_] := (
    atlasA = atlasData[[1]];
    polyAtlasA = atlasData[[2]];
    pairsAtlasA = atlasData[[3]];
    atlasB = atlasData[[4]];
    polyAtlasB = atlasData[[5]];
    pairsAtlasB = atlasData[[6]];
    
    outputData = {};

    For[i = 1, i < Length[pairs]+1, i++,
        If[Mod[i,1000] == 0, Print[i, " / ", Length[pairs], " pairs processed."], Null];
        
        pair = pairs[[i]];

        poly1 = ToExpression[poly[[pair[[1]]]][[1]]];
        poly2 = ToExpression[poly[[pair[[2]]]][[1]]];

        (* Find the possible a and b values of intersections as their indices in the atlas. *)
        aValues = {};
        For[ j = 1, j < Length[pairsAtlasA]+1, j++,
            If[ MemberQ[ pairsAtlasA[[j]], pair ], AppendTo[aValues, j], None ]
        ];
        
        bValues = {};
        For[ j = 1, j < Length[pairsAtlasB]+1, j++,
            If[ MemberQ[ pairsAtlasB[[j]], pair ], AppendTo[bValues, j], None ]
        ];

        intersections = {};
        (* Check which pairs (a value, b value) could be intersections. *)
        For[ j = 1, j < Length[aValues]+1, j++,
            aIndex = aValues[[j]];
            aVal = atlasA[[aIndex]];
            aAcc = Floor[ Accuracy[aVal] ];

            aValLowerBound = Round[aVal, 10^(-aAcc)] - 10^(-aAcc);
            aValUpperBound = Round[aVal, 10^(-aAcc)] + 10^(-aAcc);

            For[ k = 1, k < Length[bValues]+1, k++,
                bIndex = bValues[[k]];
                bVal = atlasB[[bIndex]];
                bAcc = Floor[ Accuracy[bVal] ];

                bValLowerBound = Round[bVal, 10^(-bAcc)] - 10^(-bAcc);
                bValUpperBound = Round[bVal, 10^(-bAcc)] + 10^(-bAcc);

                (* Here there should be two cases - either the result of Reduce is False, in which case there are no solutions,
                or Reduce gives the value of exactly 1 intersection. *)
                red = Reduce[And[poly1 == 0, poly2 == 0, aValLowerBound < a < aValUpperBound, bValLowerBound < b < bValUpperBound] ];
                
                If[ red === False,
                    Null,
                    If[ Or[MatchQ[red, a == _ && b == _], MatchQ[red, b == _ && a == _] ], (* If True, Reduce gave the solution we expected. *)
                    AppendTo[intersections, {aIndex, bIndex}],
                    Print["Warning - found unexpected answer when reducing at pair", i, ". red = ", red];
                    ];
                ];
            ];
        ];

        If[Length[intersections] > 0,
            AppendTo[outputData, {pair, intersections}],
        Null];
    ];

    outputData
)

(* sT orbit type *)
Print["Processing sT orbit type..."]
sTIntersections = FindIntersections[pairsT, polysT, sTAtlasData]
Export["data/sTIntersections.txt", sTIntersections]

(* gT orbit type *)
Print["Processing gT orbit type..."]
gTIntersections = FindIntersections[pairgT, polygT, gTAtlasData]
Export["data/gTIntersections.txt", gTIntersections]

(* gP orbit type *)
Print["Processing gP orbit type..."]
gPIntersections = FindIntersections[pairgP, polygP, gPAtlasData]
Export["data/gPIntersections.txt", gPIntersections]

(* sC orbit type *)
Print["Processing sC orbit type..."]
sCIntersections = FindIntersections[pairsC, polysC, sCAtlasData]
Export["data/sCIntersections.txt", sCIntersections]

(* gC orbit type *)
Print["Processing gC orbit type..."]
gCIntersections = FindIntersections[pairgC, polygC, gCAtlasData]
Export["data/gCIntersections.txt", gCIntersections]

(* sD orbit type *)
Print["Processing sD orbit type..."]
sDIntersections = FindIntersections[pairsD, polysD, sDAtlasData]
Export["data/sDIntersections.txt", sDIntersections]

(* gD orbit type *)
Print["Processing gD orbit type..."]
gDIntersections = FindIntersections[pairgD, polygD, gDAtlasData]
Export["data/gDIntersections.txt", gDIntersections]

Print["Finished!"]

Pause[20000]