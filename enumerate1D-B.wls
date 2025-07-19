#!/usr/bin/env wolframscript

Off[General::stop]

Print["Importing..."]

(* Import polynomial data from enumerate1D-A.py *)

polytT = Import["data/tTConfPoly.txt", "Table"]
polyrT = Import["data/rTConfPoly.txt", "Table"]
polyrP = Import["data/rPConfPoly.txt", "Table"]
polytO = Import["data/tOConfPoly.txt", "Table"]
polytC = Import["data/tCConfPoly.txt", "Table"]
polyrC = Import["data/rCConfPoly.txt", "Table"]
polytI = Import["data/tIConfPoly.txt", "Table"]
polytD = Import["data/tDConfPoly.txt", "Table"]
polyrD = Import["data/rDConfPoly.txt", "Table"]

Print["Factoring small orbit types (tT,rT,rP,tO,tC,rC)..."]

(* tT orbit type *)
tTfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
tTfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polytT]+1, i++,
	poly = ToExpression[polytT[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[tTfactors, factor],
                index = Position[tTfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[tTfactors];
                AppendTo[tTfactors, factor] (* False *)];
            AppendTo[tTfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/tTConfFactors.txt", tTfactors, "Table"]
Export["data/tTConfFactorIndices.txt", tTfactorIndices, "Table"]

(* rT orbit type *)
rTfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
rTfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polyrT]+1, i++,
	poly = ToExpression[polyrT[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[rTfactors, factor],
                index = Position[rTfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[rTfactors];
                AppendTo[rTfactors, factor] (* False *)];
            AppendTo[rTfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/rTConfFactors.txt", rTfactors, "Table"]
Export["data/rTConfFactorIndices.txt", rTfactorIndices, "Table"]

(* rP orbit type *)
rPfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
rPfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polyrP]+1, i++,
	poly = ToExpression[polyrP[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[rPfactors, factor],
                index = Position[rPfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[rPfactors];
                AppendTo[rPfactors, factor] (* False *)];
            AppendTo[rPfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/rPConfFactors.txt", rPfactors, "Table"]
Export["data/rPConfFactorIndices.txt", rPfactorIndices, "Table"]

(* tO orbit type *)
tOfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
tOfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polytO]+1, i++,
	poly = ToExpression[polytO[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[tOfactors, factor],
                index = Position[tOfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[tOfactors];
                AppendTo[tOfactors, factor] (* False *)];
            AppendTo[tOfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/tOConfFactors.txt", tOfactors, "Table"]
Export["data/tOConfFactorIndices.txt", tOfactorIndices, "Table"]

(* tC orbit type *)
tCfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
tCfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polytC]+1, i++,
	poly = ToExpression[polytC[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[tCfactors, factor],
                index = Position[tCfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[tCfactors];
                AppendTo[tCfactors, factor] (* False *)];
            AppendTo[tCfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/tCConfFactors.txt", tCfactors, "Table"]
Export["data/tCConfFactorIndices.txt", tCfactorIndices, "Table"]

(* rC orbit type *)
rCfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
rCfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polyrC]+1, i++,
	poly = ToExpression[polyrC[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[rCfactors, factor],
                index = Position[rCfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[rCfactors];
                AppendTo[rCfactors, factor] (* False *)];
            AppendTo[rCfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/rCConfFactors.txt", rCfactors, "Table"]
Export["data/rCConfFactorIndices.txt", rCfactorIndices, "Table"]

Print["Factoring tI orbit type..."]

(* tI orbit type *)
tIfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
tIfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polytI]+1, i++,
	poly = ToExpression[polytI[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[tIfactors, factor],
                index = Position[tIfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[tIfactors];
                AppendTo[tIfactors, factor] (* False *)];
            AppendTo[tIfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/tIConfFactors.txt", tIfactors, "Table"]
Export["data/tIConfFactorIndices.txt", tIfactorIndices, "Table"]

Print["Factoring tD orbit type..."]

(* tD orbit type *)
tDfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
tDfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polytD]+1, i++,
	poly = ToExpression[polytD[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[tDfactors, factor],
                index = Position[tDfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[tDfactors];
                AppendTo[tDfactors, factor] (* False *)];
            AppendTo[tDfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/tDConfFactors.txt", tDfactors, "Table"]
Export["data/tDConfFactorIndices.txt", tDfactorIndices, "Table"]

Print["Factoring rD orbit type..."]

(* rD orbit type *)
rDfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
rDfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polyrD]+1, i++,
	poly = ToExpression[polyrD[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[rDfactors, factor],
                index = Position[rDfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[rDfactors];
                AppendTo[rDfactors, factor] (* False *)];
            AppendTo[rDfactorIndices, {i-1, index}];
        ]
        ];
    ]

Export["data/rDConfFactors.txt", rDfactors, "Table"]
Export["data/rDConfFactorIndices.txt", rDfactorIndices, "Table"]

Print["Finished!"]

Pause[20]