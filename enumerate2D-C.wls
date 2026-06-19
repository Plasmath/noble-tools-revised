(* enumerate2D-C.wls *)
(* Code that determines the atlases of roots for each orbit type. Each root is identified by an approximate value (an upper & lower bound)
   which does not overlap with any other root in the atlas, as well as a minimal polynomial which allows for the root's value to be 
   determined to any level of precision. *)

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

(* Determines if two roots in a set of roots are close enough that we can't determine for sure if they're distinct. *)
HasCloseRoots[roots_] := (
    isTrue = False;
    For[ ii = 1, ii < Length[roots]+1, ii++,
        rootI = roots[[ii]]; (* first root for comparison *)
        For[ jj = ii+1, jj < Length[roots]+1, jj++,
            rootJ = roots[[jj]]; (* second root for comparison *)
            If[rootI == rootJ, isTrue = True, Null]; (* == evaluates to true if their upper and lower bounds overlap  *)
        ]
    ];
    isTrue
)

ProcessPairs[pairs_, poly_] := (
    atlasA = {};
    polyAtlasA = {}; (* List of polynomials in the parameter atlas, used to determine if two roots are actually the same value or just very close. *)
    pairsAtlasA = {}; (* List of pairs which contain a certain value of a at a given location in the atlas. *)
    atlasB = {};
    polyAtlasB = {}; (* Same for b *)
    pairsAtlasB = {}; (* Same for b *)

    For[i = 1, i < Length[pairs]+1, i++,
        If[Mod[i,1000] == 0, Print[i, " / ", Length[pairs], " pairs processed."], Null]; (* Periodic progress updates *)
        
        ind1 = pairs[[i]][[1]];
        ind2 = pairs[[i]][[2]];

        s = ToExpression[poly[[ind1]][[1]]]; (* First critical curve *)
        t = ToExpression[poly[[ind2]][[1]]]; (* Second critical curve *)

        (* We factorize these and find roots individually as otherwise the root-finding algorithm can have problems with certain large polynomials. *)
        resA = FactorList[Resultant[s,t,a], Extension -> {Sqrt[2],Sqrt[5]} ];
        resB = FactorList[Resultant[s,t,b], Extension -> {Sqrt[2],Sqrt[5]} ];

        rootsA = {};
        For[j = 2, j < Length[resB]+1, j++, (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
            (* Add roots of factor to full list of roots. *)  
            factor = resB[[j]][[1]] == 0; (* Factors are formatted in the form {polynomial, multiplicity}. Including only the first part of the
                                            factor stops us from needing to consider roots with multiplicity, which is very good for the creation 
					    of an atlas of roots as it ensures that each root of this polynomial we are adding to the atlas is a
					    unique root of the polynomial. *)
            If[factor == False, Continue[], Null]; (*Again, we do not consider constant factors.*)

            newRoots = a /. {ToRules[ NRoots[factor, a, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

            (* If there exist two roots of a factor that are very close to each other (i.e. upper and lower bounds overlap),
	       we must increase the precision until they are distinct and no longer overlap. *)
            precisionAmount = 1;
            While[HasCloseRoots[newRoots],
                precisionAmount++;
                newRoots = a /. {ToRules[ NRoots[factor, a, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount, Method -> "JenkinsTraub"] ]};
            ];

            rootsA = Join[rootsA, Table[{root, factor}, {root, newRoots}] ]; (* Add new roots to our existing list *)
        ];
        rootsA = Cases[rootsA, {root_, _} /; And[ PossibleZeroQ[ Im[root] ], Re[root] =!= 0, Negative[ Re[root] ] =!= True ] ]; (* Filter out roots which we have determined must be complex or negative. *)

        (* Repeat for b *)
        rootsB = {};
        For[j = 2, j < Length[resA]+1, j++,
            factor = resA[[j]][[1]] == 0;
            If[factor == False, Continue[], Null];

            newRoots = b /. {ToRules[ NRoots[factor, b, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

            precisionAmount = 1;
            While[HasCloseRoots[newRoots],
                precisionAmount++;
                newRoots = b /. {ToRules[ NRoots[factor, b, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount, Method -> "JenkinsTraub"] ]};
            ];

            rootsB = Join[rootsB, Table[{root, factor}, {root, newRoots}] ];
        ];
        rootsB = Cases[rootsB, {root_, _} /; And[ PossibleZeroQ[ Im[root] ], Re[root] =!= 0, Negative[ Re[root] ] =!= True ] ]; (* Filter out roots which we have determined must be complex or negative. *)

        (* If either rootsA or rootsB has no positive roots, then the two curves cannot have an intersection in the area we are considering (pairs of positive reals). *)
        If[Or[Length[rootsA] == 0, Length[rootsB] == 0], Continue[], Null ];

        (* We now have a set of roots which give approximate coordinates for values of the a and b parameters.
           From here, we incorporate these roots into the two atlases. *)
        For[j = 1, j < Length[rootsA]+1, j++,
            
            root = rootsA[[j]][[1]];
            polyRoot = rootsA[[j]][[2]];
            isNewRoot = True; (* We assume it is a new root unless proven otherwise. *)

            (* First, we check to see if any roots currently in the atlas overlap with the approximation of our current root.
               If not, we may add it to the atlas as it cannot be the same root as an existing one which is there already. *)
            For[k = 1, k < Length[atlasA]+1, k++,
                atlasRoot = atlasA[[k]];

                If[ atlasRoot == root, (*Upper and lower bounds overlap, so we must check if they are the same root or not. *)
                    polyAtlasRoot = polyAtlasA[[k]];

                    (*Quickly, we check if the polynomials of the two candidate roots are equal. If this is the case, they must be the same root.*)
                    If[polyAtlasRoot == polyRoot,
                        isNewRoot = False;
                        AppendTo[pairsAtlasA[[k]], {ind1, ind2}]; (* We add the pair to the relevant list whenever we find out a given root is already in the root atlas. *)
                        Break[],
                    Null];
                    (*If this fails, we check if the gcd of the two polynomials contains this root. If so, they must be the same. If not, they are different.*)
                    gcd = PolynomialGCD[ First[polyAtlasRoot], First[polyRoot], Extension -> Automatic ] == 0;

                    If[gcd == False, Continue[], Null]; (* If the gcd is a constant then the polynomials share no common factors and this may be a new root. *)

                    If[Or[gcd == polyRoot, gcd == polyAtlasRoot], (* If the gcd is either of the two source polynomials, the roots of one must be a subset of the other and therefore we are converging on the same root. *)
                        isNewRoot = False;
                        AppendTo[pairsAtlasA[[k]], {ind1, ind2}];
                        Break[],
                    Null];

                    gcdRoots = a /. {ToRules[ NRoots[gcd, a, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

                    precisionAmount = 1;
                    While[HasCloseRoots[gcdRoots],
                        precisionAmount++;
                        gcdRoots = a /. {ToRules[ NRoots[gcd, a, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount, Method -> "JenkinsTraub"] ]};
                    ];
                    gcdRoots = Cases[gcdRoots, root_ /; And[ PossibleZeroQ[ Im[root] ], Re[root] =!= 0, Negative[ Re[root] ] =!= True ] ]; (* Filter out roots which we have determined must be complex or negative. *)

                    (* Checking to see if the gcd has the root we're looking for... *)
                    For[l = 1, l < Length[gcdRoots]+1, l++,
                        gcdRoot = gcdRoots[[l]];
                        If[And[gcdRoot == atlasRoot, gcdRoot == root], (* The gcd also contains this root, so all three roots must be the same. *)
                            isNewRoot = False;
                            Break[],
                        Null ];
                    ];

                    (* Break if gcd shares this root with the other polynomials.
                       If not, the two roots are different and we must increase
                       the precision of the two approximations so the bounds do
                       not overlap. *)
                    If[Not[isNewRoot],
                        AppendTo[pairsAtlasA[[k]], {ind1, ind2}];
                        Break[], 
                    Null]; 
                    
                    While[atlasRoot == root, (* Increase precision. *)
                        precisionAmount++;
                        newAtlasRoot = FindRoot[polyAtlasRoot, {a,atlasRoot}, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount];
                        newRoot = FindRoot[polyRoot, {a,root}, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount];

                        (* Check if the new approximation is still approximating the same root. If so, set the better approximation as the value of the root. *)
			(* The program will stall if it cannot find an accurate enough approximation, so by it successfully running we know it succeeded for every root. *)
                        If[newAtlasRoot == atlasRoot,
                            atlasRoot = newAtlasRoot,
                            Print["Warning: struggling to obtain accurate approximation for ", polyAtlasRoot]
                        ];
                        If[newRoot == root,
                            root = newRoot,
                            Print["Warning: struggling to obtain accurate approximation for ", polyRoot]
                        ];
                    ];

                    atlasA[[k]] = atlasRoot;
                Null];
            ];

            If[isNewRoot, (* Add to atlas if it has passed all tests. *)
                AppendTo[atlasA, root];
                AppendTo[polyAtlasA, polyRoot];
                AppendTo[pairsAtlasA, {{ind1, ind2}}],
            Null ];
        ];

        (* Same for b *)
        For[j = 1, j < Length[rootsB]+1, j++,
            
            root = rootsB[[j]][[1]];
            polyRoot = rootsB[[j]][[2]];
            isNewRoot = True;

            For[k = 1, k < Length[atlasB]+1, k++,
                atlasRoot = atlasB[[k]];

                If[ atlasRoot == root,
                    polyAtlasRoot = polyAtlasB[[k]];

                    If[polyAtlasRoot == polyRoot,
                        isNewRoot = False;
                        AppendTo[pairsAtlasB[[k]], {ind1, ind2}];
                        Break[],
                    Null];
                    
                    gcd = PolynomialGCD[ First[polyAtlasRoot], First[polyRoot], Extension -> Automatic ] == 0;

                    If[gcd == False, Continue[], Null];

                    If[Or[gcd == polyRoot, gcd == polyAtlasRoot],
                        isNewRoot = False;
                        AppendTo[pairsAtlasB[[k]], {ind1, ind2}];
                        Break[],
                    Null];

                    gcdRoots = b /. {ToRules[ NRoots[gcd, b, PrecisionGoal -> 5, MaxIterations -> 800, Method -> "JenkinsTraub"] ]};

                    precisionAmount = 1;
                    While[HasCloseRoots[gcdRoots],
                        precisionAmount++;
                        gcdRoots = b /. {ToRules[ NRoots[gcd, b, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount, Method -> "JenkinsTraub"] ]};
                    ];
                    gcdRoots = Cases[gcdRoots, root_ /; And[ PossibleZeroQ[ Im[root] ], Re[root] =!= 0, Negative[ Re[root] ] =!= True ] ];

                    For[l = 1, l < Length[gcdRoots]+1, l++,
                        gcdRoot = gcdRoots[[l]];
                        If[And[gcdRoot == atlasRoot, gcdRoot == root],
                            isNewRoot = False;
                            Break[],
                        Null ];
                    ];

                    If[Not[isNewRoot],
                        AppendTo[pairsAtlasB[[k]], {ind1, ind2}];
                        Break[], 
                    Null]; 
                    
                    While[atlasRoot == root,
                        precisionAmount++;
                        newAtlasRoot = FindRoot[polyAtlasRoot, {b,atlasRoot}, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount];
                        newRoot = FindRoot[polyRoot, {b,root}, PrecisionGoal -> 5*precisionAmount, MaxIterations -> 800*precisionAmount];

                        If[newAtlasRoot == atlasRoot,
                            atlasRoot = newAtlasRoot,
                            Print["Warning: struggling to obtain accurate approximation for ", polyAtlasRoot]
                        ];
                        If[newRoot == root,
                            root = newRoot,
                            Print["Warning: struggling to obtain accurate approximation for ", polyRoot]
                        ];
                    ];

                    atlasB[[k]] = atlasRoot;
                Null];
            ];

            If[isNewRoot,
                AppendTo[atlasB, root];
                AppendTo[polyAtlasB, polyRoot];
                AppendTo[pairsAtlasB, {{ind1, ind2}}],
            Null ];
        ];
    ];

    output = {atlasA, polyAtlasA, pairsAtlasA, atlasB, polyAtlasB, pairsAtlasB}; output
)

(* sT orbit type *)
Print["Processing sT orbit type..."]
sTData = ProcessPairs[pairsT, polysT]
Export["data/sTAtlas.txt", sTData]

(* gT orbit type *)
Print["Processing gT orbit type..."]
gTData = ProcessPairs[pairgT, polygT]
Export["data/gTAtlas.txt", gTData]

(* gP orbit type *)
Print["Processing gP orbit type..."]
gPData = ProcessPairs[pairgP, polygP]
Export["data/gPAtlas.txt", gPData]

(* sC orbit type *)
Print["Processing sC orbit type..."]
sCData = ProcessPairs[pairsC, polysC]
Export["data/sCAtlas.txt", sCData]

(* gC orbit type *)
Print["Processing gC orbit type..."]
gCData = ProcessPairs[pairgC, polygC]
Export["data/gCAtlas.txt", gCData]

(* sD orbit type *)
Print["Processing sD orbit type..."]
sDData = ProcessPairs[pairsD, polysD]
Export["data/sDAtlas.txt", sDData]

(* gD orbit type *)
Print["Processing gD orbit type..."]
gDData = ProcessPairs[pairgD, polygD]
Export["data/gDAtlas.txt", gDData]

Print["Finished!"]

Pause[20]