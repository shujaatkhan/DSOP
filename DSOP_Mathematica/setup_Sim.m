(* ::Package:: *)

(* ::Section:: *)
(*Simulation parameters*)


SimPeople=10000;                                                                   (* Number of people to simulate *)
SimPeriods=60-25;                                                               (* Length of life in simulation (from 25 to 60 years old) *)
(* NOTE: SAK 08/13/2016
The number of simulation periods should be 36, because 25 to 60 years (inclusive) is 36 years.*)

InitialWYRatio = {0.17, .5, .83};                              (* Initial wealth-income ratio (from Cagetti (2003) p .13) *)
(*
Subscript[\[Sigma], WY] = 1.784;  
Subscript[\[Mu], WY] = -2.794;
InitialWYRatio = Exp[Log[DiscreteMeanOneLogNormal[Subscript[\[Sigma], WY],7][[1]]]+Subscript[\[Mu], WY]+1/2*Subscript[\[Sigma], WY]^2];
(* JYao: richer intial WY ratio (from Gourinchas and Parker 2002, Table II). 
We assume initial WY ratio follows a lognormal distribution with mean Subscript[\[Mu], WY] and 
standard deviation Subscript[\[Sigma], WY]. We then discretize the distribution into equiprobable 
5 values *)
*)
If[NoMoM==True,
OptionConsFnInSim={"cInterpNoMoM"};,
(* This is the default option. Method of moderation is NOT used. Endogeneous grid points are used.*)
OptionConsFnInSim={"cVecMoM", "\[Chi]Interp"};]
(* This is the default option
, on which consumption functions will be used in Simulate. *)
(* Details are in functions_ConNVal.m. *)

(*
(* begin{CDCPrivate} *)
OptionConsFnInSim={"cListable", "\[Chi]Listable"};
(* With the TighterUpBd refinement, compute consumptions vectorwise is impossible, 
and we have to do one by one. Hence we repy on the two listable functions.
Details are in functions_ConNValTighterUpBd.m. *)
(* end{CDCPrivate} *)
*)


(* ::Section:: *)
(*Drawing shocks for simulation*)


(* "ConstructShockDistribution" constructs shock distribution (JYao: varying for each time of simulation during working life) *)  
\[Theta]AgeSim[t_] := Block[{\[Sigma]\[Theta]Age,\[Theta]Age,\[Theta]Sim,\[Theta]SimProb},
\[Sigma]\[Theta]Age= \[Sigma]\[Theta][[t]];
UnemployedPeople   = SimPeople \[Mho];
{\[Theta]Sim, \[Theta]SimProb} = DiscreteMeanOneLogNormal[\[Sigma]\[Theta]Age,SimPeople - UnemployedPeople];
\[Theta]Sim=\[Theta]Sim/(1 - \[Mho]);
\[Theta]Sim   = Join[\[Theta]Sim, Table[0,{UnemployedPeople}]];
\[Theta]SimProb = \[Theta]SimProb (1 - \[Mho]) ;
\[Theta]SimProb = Join[\[Theta]SimProb, Table[\[Mho]/UnemployedPeople,{UnemployedPeople}]];
\[Theta]Age = {\[Theta]Sim, \[Theta]SimProb};
Return[\[Theta]Age];
];

\[CapitalPsi]AgeSim[t_] := Block[{\[Sigma]\[CapitalPsi]Age,\[CapitalPsi]Age},
\[Sigma]\[CapitalPsi]Age = \[Sigma]\[CapitalPsi][[t]];
\[CapitalPsi]Age = DiscreteMeanOneLogNormal[\[Sigma]\[CapitalPsi]Age,SimPeople];
Return[\[CapitalPsi]Age];
];    

(* "ConstructSimShocks" constructs simulation shocks for each period and each individual and assigns initial wealth holdings *)
ConstructSimShocks := Block[{},
If[Sdfix==True, (* If the st. dev. is fixed, then we can speed up the code by drawing from the distribution just once *)
If[SimPeriods>40,Print["Simulating retirement years as well. Write the code for retirement year draws."];Abort[]];
UnemployedPeople   = SimPeople \[Mho];
Sdfix\[Theta] = \[Sigma]\[Theta][[1]];
Sdfix\[CapitalPsi] = \[Sigma]\[CapitalPsi][[1]];
{\[Theta]SimFix, \[Theta]SimProbFix} = DiscreteMeanOneLogNormal[Sdfix\[Theta],SimPeople - UnemployedPeople];
\[Theta]SimFix=\[Theta]SimFix/(1 - \[Mho]); \[Theta]SimFix = Join[\[Theta]SimFix, Table[0,{UnemployedPeople}]];
\[Theta]SimProbFix = \[Theta]SimProbFix (1 - \[Mho]) ; \[Theta]SimProbFix = Join[\[Theta]SimProbFix, Table[\[Mho]/UnemployedPeople,{UnemployedPeople}]];
{\[CapitalPsi]SimFix, \[CapitalPsi]SimProbFix} = DiscreteMeanOneLogNormal[Sdfix\[CapitalPsi],SimPeople];
\[Theta]SimMat=Table[RandomSample[\[Theta]SimFix],{t,1,SimPeriods}];
\[CapitalPsi]SimMat=Table[RandomSample[\[CapitalPsi]SimFix],{t,1,SimPeriods}];
,
\[Theta]SimMat=Table[RandomSample[\[Theta]AgeSim[t][[1]]],{t,1,SimPeriods}];
\[CapitalPsi]SimMat=Table[RandomSample[\[CapitalPsi]AgeSim[t][[1]]],{t,1,SimPeriods}];
];
w0Sim=RandomChoice[InitialWYRatio,SimPeople]; ];               (* Initial wealth to income ratio *)


(* ::Section:: *)
(*Simulated medians*)


(* ::Subsection:: *)
(*Using CPU*)


(* Computing simulated medians *)
Simulate := Block[{},
Do[FuncToDefine\[ScriptC]VecN\[Chi]Vec[OptionConsFnInSim], {1}];
ClearAll[mtSimMatCPU, ctSimMatCPU, atSimMatCPU, wtSimMatCPU];

(* NOTE: The simulation needs to happen for age 25 to 60 (inclusive).
         \[Theta]SimMat, however, starts at age 26. To simulate age 25, we
         assume that everyone gets a \[Theta] shock equal to 1.
*)
mtSimMat= ctSimMat=atSimMat={};
(* Simulating age 25 *)
wtSimMat={w0Sim};
AppendTo[mtSimMat, wtSimMat[[-1]]+ConstantArray[1,SimPeople]]; (* \[Theta]=1 at a=25*)
AppendTo[ctSimMat, \[ScriptC]Vec[mtSimMat[[-1]], PeriodsToSolve]]; 
AppendTo[atSimMat, (mtSimMat[[-1]]-ctSimMat[[-1]])];

Table[
AppendTo[wtSimMat,\[DoubleStruckCapitalR]/(\[ScriptCapitalG]Vect[[t]]\[CapitalPsi]SimMat[[t]])*atSimMat[[-1]]];
AppendTo[mtSimMat, wtSimMat[[-1]] + \[Theta]SimMat[[t]]];
AppendTo[ctSimMat, \[ScriptC]Vec[mtSimMat[[-1]], PeriodsToSolve-t]]; 
AppendTo[atSimMat, (mtSimMat[[-1]]-ctSimMat[[-1]])]; 
,{t,1,SimPeriods}]; (* We are simulating from age 25 to 60 (inclusive)*)
wtSimMat=Drop[wtSimMat,1]; (* Dropping the initial wealth to income, i.e. age 25*)

SimMedian= Table[Median[Flatten[Table[wtSimMat[[5 (t - 1) + s]], {s, 1, 5}], 1]], {t, 1, 7}];

(* The following three lines are added for comparing CPU and GPU results. *)
If[CPUVsGPUTestOption==True
,
mtSimMatCPU=mtSimMat;
ctSimMatCPU=ctSimMat;
atSimMatCPU=atSimMat;
wtSimMatCPU=wtSimMat;
];
ClearAll[mtSimMat, ctSimMat, atSimMat, wtSimMat];
];


(* ::Subsection:: *)
(*Using GPU*)


CPUVsGPUTestOption=True;


(* MNW: Computes simulated medians using the GPU routine below. *)
SimulateGPU := Block[{},
Do[FuncToDefine\[ScriptC]VecN\[Chi]Vec[OptionConsFnInSim], {1}];
ClearAll[mtSimMatGPU, ctSimMatGPU, atSimMatGPU, wtSimMatGPU];

{mtSimMat, ctSimMat, atSimMat, wtSimMat} = GPUsimFunc[w0Sim];
SimMedian= Table[Median[Flatten[Table[wtSimMat[[5 (t - 1) + s]], {s, 1, 5}], 1]], {t, 1, 7}];

(* The following three lines are added for comparing CPU and GPU results. *)
If[CPUVsGPUTestOption==True
,
mtSimMatGPU=mtSimMat;
ctSimMatGPU=ctSimMat;
atSimMatGPU=atSimMat;
wtSimMatGPU=wtSimMat;
];
ClearAll[mtSimMat, ctSimMat, atSimMat, wtSimMat];

];
