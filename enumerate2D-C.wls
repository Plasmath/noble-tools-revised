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

Print["Processing pair data..."]

(* Function for pair-processing, used for all orbit types. Returns data to be exported. *)
processPairs[pairs_, poly_] := (
data = {};
For[i = 1, i < Length[pairs]+1, i++,
    If[Mod[i,2000] == 0, Print[i, " / ", Length[pairs], " pairs processed."], Null];
    
    ind1 = pairs[[i]][[1]];
    ind2 = pairs[[i]][[2]];

    s = ToExpression[poly[[ind1]][[1]]];
    t = ToExpression[poly[[ind2]][[1]]];

    (* We factorize these and find roots individually as otherwise the root-finding algorithm can have problems with certain large polynomials. *)
    resA = FactorList[Resultant[s,t,a], Extension -> Sqrt[5] ];
    resB = FactorList[Resultant[s,t,b], Extension -> Sqrt[5] ];

    rootsA = {};
    For[j = 2, j < Length[resB]+1, j++, (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
        (* Add roots of factor to full list of roots. *)   
        factor = resB[[j]][[1]] == 0;
        If[factor == False, Continue[], Null]; (*Again, we do not consider constant factors.*)
        
        newRoots = a /. {ToRules[ NRoots[factor, a, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

        rootsA = Join[rootsA, newRoots];
    ];

    rootsB = {};
    For[j = 2, j < Length[resA]+1, j++, (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
        (* Add roots of factor to full list of roots. *)   
        factor = resA[[j]][[1]] == 0;
        If[factor == False, Continue[], Null]; (*Again, we do not consider constant factors.*)
        
        newRoots = b /. {ToRules[ NRoots[factor, b, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

        rootsB = Join[rootsB, newRoots];
    ];

    (* When either resultant is a constant, one obtains a boolean value rather than an equation.
    As these trivially cannot have roots, we may skip this pair.*)
    
    For[j = 1, j < Length[rootsA]+1, j++,
        root = rootsA[[j]];

        If[ PossibleZeroQ[ Im[root] ], (*Check if root precision includes the possibility of this being a real root*)
            root = Re[root];
            If[ Positive[root], (*Check if root is positive, given that it should be real*)
                AppendTo[ data, {ind1, ind2, a, ToString[root, InputForm]} ], (*Indicate existence of root for a*)
                Null,
                AppendTo[ data, {ind1, ind2, a, ToString[root, InputForm]} ] (*Still do this in the case of an inconclusive result.*)
            ],
            Null
        ];
    ];

    For[j = 1, j < Length[rootsB]+1, j++, (*Same process for b*)
        root = rootsB[[j]];
        If[ PossibleZeroQ[ Im[root] ],
            root = Re[root];
            If[ Positive[root],
                AppendTo[ data, {ind1, ind2, b, ToString[root, InputForm]} ],
                Null,
                AppendTo[ data, {ind1, ind2, b, ToString[root, InputForm]} ]
            ],
            Null
        ];
    ];
]; data)

(* sT orbit type *)
Print["Processing sT orbit type..."]
sTData = processPairs[pairsT, polysT]
Export["data/sTRootData.txt",sTData,"Table"]

(* gT orbit type *)
Print["Processing gT orbit type..."]
gTData = processPairs[pairgT, polygT]
Export["data/gTRootData.txt",gTData,"Table"]

(* gP orbit type *)
Print["Processing gP orbit type..."]
gPData = processPairs[pairgP, polygP]
Export["data/gPRootData.txt",gPData,"Table"]

(* sC orbit type *)
Print["Processing sC orbit type..."]
sCData = processPairs[pairsC, polysC]
Export["data/sCRootData.txt",sCData,"Table"]

(* gC orbit type *)
Print["Processing gC orbit type..."]
gCData = processPairs[pairgC, polygC]
Export["data/gCRootData.txt",gCData,"Table"]

(* sD orbit type *)
Print["Processing sD orbit type..."]
sDData = processPairs[pairsD, polysD]
Export["data/sDRootData.txt",sDData,"Table"]

(* gD orbit type *)
Print["Processing gD orbit type..."]
gDData = processPairs[pairgD, polygD]
Export["data/gDRootData.txt",gDData,"Table"]

Print["Finished!"]

Pause[20]