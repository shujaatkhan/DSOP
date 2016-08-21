(* ::Package:: *)

ClearAll[FuncToDefine\[GothicV]Func];
FuncToDefine\[GothicV]Func[\[Rho]_?NumericQ]:=Block[{},

Clear[\[GothicV],\[GothicV]a,\[GothicV]aa,\[GothicC],\[GothicC]a];
(* The idea here is that with consumption function, MPC function and value function at t+1,
we could derive the GothicV function and its first two derivatives at t.  *)

(* Note the difference between this section and the previous 2period: here we have 
shocks to both permanent and transitory income. *)

\[GothicV][at_?NumericQ,PeriodsUntilT_] := Block[{mtp1, \[GothicV]Val},
  If[PeriodsUntilT == 0,(*then*)Return[0.]];
  (*else*)
  \[GothicV]Val=\[Beta]Life[[PeriodsUntilT+1]]*Sum[
       mtp1=at*RLife[[PeriodsUntilT+1]]/(\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]]) + \[Theta]Vals[[\[Theta]Loop]];
      Times[\[Theta]Probs[[\[Theta]Loop]]*\[CapitalPsi]Probs[[\[CapitalPsi]Loop]]
, (\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]])^(1-\[Rho])
, \[ScriptV][mtp1,PeriodsUntilT-1]
], {\[Theta]Loop,1,Length[\[Theta]Vals]}, {\[CapitalPsi]Loop,1,Length[\[CapitalPsi]Vals]}];  
  Return[\[GothicV]Val];
]; (* End Block[]*)

\[GothicV]a[at_?NumericQ,PeriodsUntilT_] := Block[{mtp1,\[GothicV]aVal},
  If[PeriodsUntilT == 0,(*then*)Return[0.]];
  (*else*)
  \[GothicV]aVal=\[Beta]Life[[PeriodsUntilT+1]]*RLife[[PeriodsUntilT+1]]*Sum[
      mtp1=at*RLife[[PeriodsUntilT+1]]/(\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]]) + \[Theta]Vals[[\[Theta]Loop]];
      Times[\[Theta]Probs[[\[Theta]Loop]]*\[CapitalPsi]Probs[[\[CapitalPsi]Loop]]
, (\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]])^(-\[Rho])
, uP[\[ScriptC][mtp1, PeriodsUntilT-1]]
], {\[Theta]Loop,1,Length[\[Theta]Vals]}, {\[CapitalPsi]Loop,1,Length[\[CapitalPsi]Vals]}];  
  Return[\[GothicV]aVal];
];(*End Block[]*)


\[GothicC][at_?NumericQ, PeriodsUntilT_] := nP[\[GothicV]a[at,PeriodsUntilT]];

SetAttributes[{\[GothicV],\[GothicV]a,\[GothicC]}, Listable]; (* Allows funcs to operate on lists *)


(* Vector-based Functions *)
ClearAll[\[DoubleStruckCapitalE]All];

\[DoubleStruckCapitalE]All[aVect_?VectorQ, PeriodsUntilT_] := Block[{aVectLength, aVectLengthZero,mtp1Vec},
aVectLength=Length[aVect];
If[PeriodsUntilT == 0
,(*then*)aVectLengthZero=Table[0., {aVectLength}]
; \[DoubleStruckCapitalE]AllThis={aVectLengthZero}];

If[PeriodsUntilT > 0, (*then*)
  mtp1Vec=Table[aVect*RLife[[PeriodsUntilT+1]]/(\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]]) + \[Theta]Vals[[\[Theta]Loop]]
               , {\[Theta]Loop,1,Length[\[Theta]Vals]}, {\[CapitalPsi]Loop,1,Length[\[CapitalPsi]Vals]}];
  (* Different values of mtp1 for different values of \[Theta], \[CapitalPsi], and aVect *)
  mtp1VecFlatten=Flatten[mtp1Vec]; (*Length[\[Theta]Vals]xLength[\[CapitalPsi]Vals]xLength[aVect]*)
  
  (* Over here we need a function that calculates c_t+1, using interpolation data from t+1 for 
     c and m, and m_t+1 
     Replace \[ScriptC]N\[Kappa]Vec with interp1 type function *)
  {\[ScriptC]VectFlatten}=Transpose[\[ScriptC]NoMoMVec[mtp1VecFlatten, PeriodsUntilT-1]];
  {\[ScriptC]Vect}=Map[Partition[Partition[#, Length[aVect]], Length[\[CapitalPsi]Vals]] &, {\[ScriptC]VectFlatten}];
  
  (* GothicVa_EG.m *)
  \[GothicV]aVal=\[Beta]Life[[PeriodsUntilT+1]]*RLife[[PeriodsUntilT+1]]*Sum[
      Times[\[Theta]Probs[[\[Theta]Loop]]*\[CapitalPsi]Probs[[\[CapitalPsi]Loop]]
      , (\[CapitalPsi]Vals[[\[CapitalPsi]Loop]]*\[ScriptCapitalG]Life[[PeriodsUntilT+1]])^(-\[Rho])
      , uP[\[ScriptC]Vect[[\[Theta]Loop, \[CapitalPsi]Loop]]]
      ], {\[Theta]Loop,1,Length[\[Theta]Vals]}, {\[CapitalPsi]Loop,1,Length[\[CapitalPsi]Vals]}];  
];

  \[DoubleStruckCapitalE]AllThis={\[GothicV]aVal};

Return[\[DoubleStruckCapitalE]AllThis];
];

]; (* End of Block *)
