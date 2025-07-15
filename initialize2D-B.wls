#!/usr/bin/env wolframscript

Off[General::stop]

Print["Importing..."]

(* Import polynomial data from initialize2D-A.wls *)

polysT = Import["data/sTConfPoly.txt", "Table"]
polygT = Import["data/gTConfPoly.txt", "Table"]
polygP = Import["data/gPConfPoly.txt", "Table"]
polysC = Import["data/sCConfPoly.txt", "Table"]
polygC = Import["data/gCConfPoly.txt", "Table"]
polysD = Import["data/sDConfPoly.txt", "Table"]
polygD = Import["data/gDConfPoly.txt", "Table"]

Print["Factoring small orbit types (sT, gT, gP, sC, gC)..."]

(* sT orbit type *)
sTfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
sTfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polysT]+1, i++,
	poly = ToExpression[polysT[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[sTfactors, factor],
                index = Position[sTfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[sTfactors];
                AppendTo[sTfactors, factor] (* False *)];
            AppendTo[sTfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/sTConfFactors.txt", sTfactors, "Table"]
Export["data/sTConfFactorIndices.txt", sTfactorIndices, "Table"]

(* gT orbit type *)
gTfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
gTfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polygT]+1, i++,
	poly = ToExpression[polygT[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[gTfactors, factor],
                index = Position[gTfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[gTfactors];
                AppendTo[gTfactors, factor] (* False *)];
            AppendTo[gTfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/gTConfFactors.txt", gTfactors, "Table"]
Export["data/gTConfFactorIndices.txt", gTfactorIndices, "Table"]

(* gP orbit type *)
gPfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
gPfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polygP]+1, i++,
	poly = ToExpression[polygP[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)

    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[gPfactors, factor],
                index = Position[gPfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[gPfactors];
                AppendTo[gPfactors, factor] (* False *)];
            AppendTo[gPfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/gPConfFactors.txt", gPfactors, "Table"]
Export["data/gPConfFactorIndices.txt", gPfactorIndices, "Table"]

(* sC orbit type *)
sCfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
sCfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polysC]+1, i++,
	poly = ToExpression[polysC[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[sCfactors, factor],
                index = Position[sCfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[sCfactors];
                AppendTo[sCfactors, factor] (* False *)];
            AppendTo[sCfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/sCConfFactors.txt", sCfactors, "Table"]
Export["data/sCConfFactorIndices.txt", sCfactorIndices, "Table"]

(* gC orbit type *)
gCfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
gCfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polygC]+1, i++,
	poly = ToExpression[polygC[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[gCfactors, factor],
                index = Position[gCfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[gCfactors];
                AppendTo[gCfactors, factor] (* False *)];
            AppendTo[gCfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/gCConfFactors.txt", gCfactors, "Table"]
Export["data/gCConfFactorIndices.txt", gCfactorIndices, "Table"]

Print["Factoring sD orbit type..."]

(* sD orbit type *)
sDfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
sDfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polysD]+1, i++,
	poly = ToExpression[polysD[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[sDfactors, factor],
                index = Position[sDfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[sDfactors];
                AppendTo[sDfactors, factor] (* False *)];
            AppendTo[sDfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/sDConfFactors.txt", sDfactors, "Table"]
Export["data/sDConfFactorIndices.txt", sDfactorIndices, "Table"]

Print["Factoring gD orbit type..."]

(* gD orbit type *)
gDfactors = {} (* List of factors under extension Q[sqrt(2),sqrt(5)] as these are the only roots that appear in configuration polynomials *)
gDfactorIndices = {} (* Indices that indicate which factors are factors of which polynomials in the volume configuration *)

For[i = 1, i < Length[polygD]+1, i++,
	poly = ToExpression[polygD[[i]][[1]]];
    polyFactors = FactorList[poly, Extension -> {Sqrt[2],Sqrt[5]}]; (* Factor under field extension *)
    
    (* We start at 2 because the first factor in the list is always a constant, which we don't consider as these do not create critical equivalence classes. *)
    For[j = 2, j < Length[polyFactors]+1, j++, 
        factor = polyFactors[[j]][[1]];
        If[NumericQ[factor], Null (* Do nothing if true *) ,
            If[MemberQ[gDfactors, factor],
                index = Position[gDfactors, factor, 1][[1]][[1]]-1 (* True *),
                index = Length[gDfactors];
                AppendTo[gDfactors, factor] (* False *)];
            AppendTo[gDfactorIndices, {i-1, index}];
        ]
        ];
    ]
Export["data/gDConfFactors.txt", gDfactors, "Table"]
Export["data/gDConfFactorIndices.txt", gDfactorIndices, "Table"]

Print["Finished!"]

Pause[20]