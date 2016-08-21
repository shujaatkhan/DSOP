(* ::Package:: *)

Needs["OpenCLLink`"];

(* Automatically choose the device with the most cores *)
CoreCountList = Table[Table[OpenCLInformation[p,d,"Core Count"],{d,1,Length[OpenCLInformation[p]]-4}],{p,1,Length[OpenCLInformation[]]}];
MaxCoreLocations = Position[CoreCountList,Max[CoreCountList]];
MyOpenCLPlatform = MaxCoreLocations[[1,1]];
MyOpenCLDevice = MaxCoreLocations[[1,2]];

(* Manually set platform and device if preferred *)
(*MyOpenCLPlatform = 1;
MyOpenCLDevice = 1;*) 

(* Determine whether double precision can be used *)
UseDoublePrecision = MemberQ[OpenCLInformation[MyOpenCLPlatform,MyOpenCLDevice,"Extensions"],"cl_khr_fp64"] || MemberQ[OpenCLInformation[MyOpenCLPlatform,MyOpenCLDevice,"Extensions"],"cl_amd_fp64"];
TypeWordA = If[UseDoublePrecision,"double","float"]; (* for use in the kernels *)
TypeWordB = If[UseDoublePrecision,"Double","Float"]; (* for use by OpenCLFunctionLoad *)


(* MNW: This reproduces the work of Construct\[ScriptC]FuncLife, but (mostly) does it on the GPU *)
Construct\[ScriptC]FuncLifeGPU[\[Rho]_,\[Bet]_]:= Block[{\[Mu]GridLifeOut,CoeffsLifeOut},

(* Generate life cycle lists to be passed to the GPU *)
\[Beta]LifeYoungToOld = Table[\[Bet]*\[Beta]hat[[i]]*\[ScriptCapitalD]Cancel[[i]], {i, 1, PeriodsToSolve}];
Do[FuncToSetupLastPeriod[], {1}];
Do[FuncToDefineUtility[\[Rho]], {1}];
Do[FuncToDefine\[GothicV]Func[\[Rho]], {1}];
Do[FuncToDefine\[CapitalLambda]Func[\[Rho]], {1}];
PeriodsSolved = 0;
Do[\[Theta]Probs = \[Theta]Mat[[PeriodsSolved+1, 2]];
   AddNewPeriodToParamLifeDates[\[Rho]];
   PeriodsSolved++;,
   {PeriodsToSolve}];
aVecLife = Table[Sort[Join[{0.03},\[GothicA]Vec]+\[GothicA]LowerBoundLife[[t+1]]],{t,1,PeriodsToSolve}];

{CoeffsLifeOut,\[Mu]GridLifeOut} = GPUkernelSolveLifeCycle[
Table[0,{(\[GothicA]GridN + 2)*4*PeriodsToSolve}] (* CoeffsLifeOut *)
, Table[0,{(\[GothicA]GridN + 1)*PeriodsToSolve}] (* \[Mu]GridLifeOut *)
, Flatten[aVecLife]
, \[Rho]
, RLife
, \[ScriptCapitalG]Life
, \[Beta]Life
, Flatten[\[Theta]ValsLife]
, Flatten[\[CapitalPsi]ValsLife]
, Flatten[\[Theta]ProbsLife]
, Flatten[\[CapitalPsi]ProbsLife]
, \[GothicA]LowerBoundLife
, \[FilledUpTriangle]\[GothicH]AccessibleLife
, \[Kappa]MinLife
, Table[0,{4*(\[GothicA]GridN+2)}] (* Current coefficients *)
, Table[0,{\[GothicA]GridN+1}] (* Current mu grid *)
, Table[0,{(\[GothicA]GridN + 1)*(NumOf\[Theta]ShockPoints+1)*NumOf\[CapitalPsi]ShockPoints}] (* TheseGothicvP *)
, Table[0,{(\[GothicA]GridN + 1)*(NumOf\[Theta]ShockPoints+1)*NumOf\[CapitalPsi]ShockPoints}] (* TheseGothicvPP *)
, Table[0,{(\[GothicA]GridN + 1)*(NumOf\[Theta]ShockPoints+1)}] (* TempGothicvP *)
, Table[0,{(\[GothicA]GridN + 1)*(NumOf\[Theta]ShockPoints+1)}] (* TempGothicvP *)
, Table[0,{\[GothicA]GridN + 1}] (* chiNow *)
, Table[0,{\[GothicA]GridN + 1}] (* chiPNow *)
]; 

\[Mu]GridsLife = Join[\[Mu]GridsLife,ArrayReshape[\[Mu]GridLifeOut,{PeriodsToSolve,(\[GothicA]GridN + 1)}]];
\[Chi]\[Mu]CoeffsLife = Join[\[Chi]\[Mu]CoeffsLife,ArrayReshape[CoeffsLifeOut,{PeriodsToSolve,(\[GothicA]GridN + 2),4}]];
];


(* This cell defines an OpenCL kernel that solves the life cycle  *)
NGridPoints = \[GothicA]GridN + 1;

LifeCycleSolveSrc = "#ifdef USING_DOUBLE_PRECISIONQ
#pragma OPENCL EXTENSION cl_khr_fp64 : enable /* Load double precision driver for Nvidia */
#pragma OPENCL EXTENSION cl_amd_fp64 : enable /* Load double precision driver for AMD */
#endif /* USING_DOUBLE_PRECISIONQ */
/* If Mma has detected that double precision can be used, these libraries are loaded */

__kernel void LifeCycleSolve_kernel(
  __global " <> TypeWordA <> " * CoeffsLifeOut /* How is the dimension of this variable determined?  Why is it double regardless of the pragma? */
, __global " <> TypeWordA <> " * muGridLifeOut /* These variables ending in Out are the only outputs; see GPUkernelSolveLifeCycle below */
, __global " <> TypeWordA <> " * aVecLife
, " <> TypeWordA <> " rho /* Only vectors need __; for some reason it crashes if rho is preceded by __ */
, __global " <> TypeWordA <> " * RLife /* Constants go into read-only memory which (presumably) is faster */
, __global " <> TypeWordA <> " * GammaLife
, __global " <> TypeWordA <> " * betaLife
, __global " <> TypeWordA <> " * ThetaVals
, __global " <> TypeWordA <> " * PsiVals
, __global " <> TypeWordA <> " * ThetaProbs
, __global " <> TypeWordA <> " * PsiProbs
, __global " <> TypeWordA <> " * BoundLife
, __global " <> TypeWordA <> " * DeltahLife
, __global " <> TypeWordA <> " * kappaMinLife
, __global " <> TypeWordA <> " * Coeffs /* Now things that the code works with */
, __global " <> TypeWordA <> " * muGrid
, __global " <> TypeWordA <> " * TheseGothicvP
, __global " <> TypeWordA <> " * TheseGothicvPP
, __global " <> TypeWordA <> " * TempGothicvP
, __global " <> TypeWordA <> " * TempGothicvPP
, __global " <> TypeWordA <> " * chiNow
, __global " <> TypeWordA <> " * chiPNow
) {

    /* Determine whether this core should work at all */
    int aN = " <> ToString[NGridPoints] <>";
    int ThetaN = " <> ToString[NumOf\[Theta]ShockPoints + 1] <> ";
    int PsiN = " <> ToString[NumOf\[CapitalPsi]ShockPoints] <> ";
    int index = get_global_id(0); /* This is the thread's ID */
    if (index >= ThetaN*PsiN*aN) return; /* We only need enough threads to handle each combo */

    /* Declare private (thread-specific) variables */
    int t; /* t=0 is first solved period (second-to-last period of life -- last period is c=m */
    int ThetaIdx;
    int PsiIdx;
    int aIdx;
    " <> TypeWordA <> " Theta;
    " <> TypeWordA <> " Psi;
    " <> TypeWordA <> " Prob;
    " <> TypeWordA <> " R;
    " <> TypeWordA <> " beta;
    " <> TypeWordA <> " Gamma;
    " <> TypeWordA <> " Bound;    /* Bound could be replaced by Boundtp1 */
    " <> TypeWordA <> " BoundNow; /* BoundNow could be replaced by Boundt */
    " <> TypeWordA <> " kappaMin;
    " <> TypeWordA <> " kappaMinNow;
    " <> TypeWordA <> " Deltah;
    " <> TypeWordA <> " DeltahNow;
    " <> TypeWordA <> " a;
    " <> TypeWordA <> " m;
    " <> TypeWordA <> " Deltam;
    " <> TypeWordA <> " mu;
    " <> TypeWordA <> " muX; /* muX is rescaled mu for the new space */
    int ii;
    int SegmentStart; /* Which segment of the interpolated function am I on */
    " <> TypeWordA <> " b0; /* Coefficients for the approximating functions */
    " <> TypeWordA <> " b1;
    " <> TypeWordA <> " b2;
    " <> TypeWordA <> " b3;
    " <> TypeWordA <> " chi;
    " <> TypeWordA <> " chiP; /* chi prime */
    " <> TypeWordA <> " Q;
    " <> TypeWordA <> " QP;
    " <> TypeWordA <> " cOptimist;
    " <> TypeWordA <> " maxBufferSaving; /* Diff between pessimist and optimist */
    " <> TypeWordA <> " cRealist;
    " <> TypeWordA <> " kappaRealist;
    " <> TypeWordA <> " Cons;
    " <> TypeWordA <> " MPC;
    " <> TypeWordA <> " GothicvP;
    " <> TypeWordA <> " GothicvPP;
    int Rem; /* Remainder is a temporary variable */
    int BaseLoc;
    " <> TypeWordA <> " ConsNow;
    " <> TypeWordA <> " dConsda;
    " <> TypeWordA <> " MPCNow;
    " <> TypeWordA <> " mNow;
    " <> TypeWordA <> " DeltamNow;
    " <> TypeWordA <> " muNow;
    " <> TypeWordA <> " ConsOpt;
    " <> TypeWordA <> " ConsPes;
    " <> TypeWordA <> " Qnow;
    " <> TypeWordA <> " QPnow;
    " <> TypeWordA <> " Span;
    " <> TypeWordA <> " chi0;
    " <> TypeWordA <> " chi1;
    " <> TypeWordA <> " chiP0;
    " <> TypeWordA <> " chiP1;
    " <> TypeWordA <> " SpanNow;

    /* Determine which asset & shock points this core refers to */
    aIdx = index/(ThetaN*PsiN);
    Rem = index - aIdx*(ThetaN*PsiN);
    ThetaIdx = Rem/PsiN;
    PsiIdx = Rem - ThetaIdx*PsiN;

    /* Begin the time loop */
    for (t = 0; t<" <> ToString[PeriodsToSolve] <>"; t++) {

        /* Update time specific variables */
        R = RLife[t+1];
        beta = betaLife[t+1];
        Gamma = GammaLife[t+1];
        Bound = BoundLife[t];
        BoundNow = BoundLife[t+1];
        Deltah = DeltahLife[t];
        DeltahNow = DeltahLife[t+1];
        kappaMin = kappaMinLife[t];
        kappaMinNow = kappaMinLife[t+1];

        /* Grab the appropriate data for this core */
        Theta = ThetaVals[ThetaIdx + t*ThetaN];
        Psi = PsiVals[PsiIdx + t*PsiN];
        Prob = ThetaProbs[ThetaIdx + t*ThetaN]*PsiProbs[PsiIdx + t*PsiN];
        a = aVecLife[aIdx + t*aN];

        /* Calculate mu next period */
        m = a*R/(Psi*Gamma) + Theta;
        Deltam = m - Bound;
        mu = log(Deltam);

        /* If next period is the last period, consume all resources */
        if (t == 0) {
            Cons = m;
            MPC = 1;
        }
        else {
            /* Get the correct coefficients for this mu */
            SegmentStart = 0;
            for (ii = 0; ii<aN; ii++) {
                if (mu > muGrid[ii]) {
                    SegmentStart += 4;
                }
            }
            b0 = Coeffs[SegmentStart+0];
            b1 = Coeffs[SegmentStart+1];
            b2 = Coeffs[SegmentStart+2];
            b3 = Coeffs[SegmentStart+3];
            if ((SegmentStart/4 > 0) && (SegmentStart/4 < aN)) {
                Span = (muGrid[SegmentStart/4] - muGrid[SegmentStart/4 - 1]);
                muX = (mu - muGrid[SegmentStart/4 - 1])/Span;
            }
            else {
                Span = 1;
                muX = mu;
            }
            barrier(CLK_GLOBAL_MEM_FENCE); /* Wait until all threads have finished their allotted tasks */

            /* Evaluate chi(mu) and chi'(mu) */
            chi = b0 + muX*(b1 + muX*(b2 + muX*(b3)));
            chiP = (b1 + muX*(2*b2 + muX*(3*b3)))/Span;

            /* Calculate Q(mu) and Q'(mu) */
            Q = 1/(1 + exp(chi));
            QP = chiP*(Q - 1)*Q;

            /* Calculate consumption and the MPC */
            cOptimist = kappaMin*(Deltam + Deltah);
            maxBufferSaving = kappaMin*(Deltah);
            cRealist = cOptimist - maxBufferSaving*Q;
            kappaRealist = kappaMin - QP*maxBufferSaving/Deltam;
            /* The rest of the corresponding line in \[ScriptC]N\[Kappa]Vec[] is unnecessary, as it is zero */

            /* Apply the liquidity constraint */
            if (cRealist >= m) {
                Cons = m - 0.000001;
                MPC = 1;
            }
            else {
                Cons = cRealist;
                MPC = kappaRealist;
            }
        }

        /* Calculate the contribution to end-of-period expected marginal value and marginal marginal value */
        TheseGothicvP[index] = Prob*powr(Psi*Gamma*Cons,-rho);/* Everybody's contributions */
        TheseGothicvPP[index] = Prob*powr(Psi*Gamma*Cons,-rho-1)*(-rho)*MPC;
        barrier(CLK_GLOBAL_MEM_FENCE); /* Other threads are about to access these, need to be sure those mem ops are done before adding */

        /* Add up the contributions across shocks for each point in aVec */
        if (PsiIdx == 0) { /* Only threads whose permanent shock index is zero do this stuff */
            BaseLoc = aIdx*(ThetaN*PsiN) + ThetaIdx*PsiN;
            TempGothicvP[aIdx*ThetaN + ThetaIdx] = 0;
            TempGothicvPP[aIdx*ThetaN + ThetaIdx] = 0;
            for (ii = 0; ii<PsiN; ii++) { 
                TempGothicvP[aIdx*ThetaN + ThetaIdx] += TheseGothicvP[BaseLoc + ii];
                TempGothicvPP[aIdx*ThetaN + ThetaIdx] += TheseGothicvPP[BaseLoc + ii];
            }
        }
        barrier(CLK_GLOBAL_MEM_FENCE); /* Make sure you're finished before you go on! */
        if ((PsiIdx + ThetaIdx) == 0) { /* Only threads where the permanent and transitory index is zero do this */
            GothicvP = 0;
            GothicvPP = 0;
            for (ii = 0; ii<ThetaN; ii++) {
                GothicvP += TempGothicvP[aIdx*ThetaN + ii];
                GothicvPP += TempGothicvPP[aIdx*ThetaN + ii];
            }
            GothicvP = GothicvP*beta*R;
            GothicvPP = GothicvPP*beta*R*R;
        }
        barrier(CLK_GLOBAL_MEM_FENCE); /* Make sure you're finished (may not need this one) */

        if ((PsiIdx + ThetaIdx) == 0) {
            /* Transform the gothicv values into consumption, MPC, and market resources */
            ConsNow = powr(GothicvP,-1/rho);
            dConsda = GothicvPP/(-rho*powr(ConsNow,-rho-1));
            MPCNow = dConsda/(dConsda + 1);
            a = aVecLife[aIdx + t*aN];
            mNow = ConsNow + a;

            /* Transform into mu and chi */
            DeltamNow = mNow - BoundNow;
            muNow = log(DeltamNow);
            muGrid[aIdx] = muNow;
            muGridLifeOut[t*aN + aIdx] = muNow; /* Report new mu grid */
            ConsOpt = kappaMinNow*(DeltamNow + DeltahNow);
            ConsPes = kappaMinNow*DeltamNow;
            Qnow = (ConsOpt - ConsNow)/(ConsOpt - ConsPes);
            QPnow = DeltamNow*(kappaMinNow - MPCNow)/(ConsOpt - ConsPes);
            chiNow[aIdx] = log((1/Qnow) - 1);
            chiPNow[aIdx] = QPnow/((Qnow - 1)*Qnow);
        }
        barrier(CLK_GLOBAL_MEM_FENCE);

        if ((PsiIdx + ThetaIdx) == 0) {
            /* Calculate the interior interval coefficients */
            if (aIdx > 0) { 
                SpanNow = muGrid[aIdx] - muGrid[aIdx-1];
                /*JunkOut[index] = SpanNow;*/
                chi0 = chiNow[aIdx - 1];
                chi1 = chiNow[aIdx];
                chiP0 = chiPNow[aIdx - 1]*SpanNow;
                chiP1 = chiPNow[aIdx]*SpanNow;
                Coeffs[(aIdx*4) + 0] = chi0;
                Coeffs[(aIdx*4) + 1] = chiP0;
                Coeffs[(aIdx*4) + 2] = 3*(chi1 - chi0) - 2*chiP0 - chiP1;
                Coeffs[(aIdx*4) + 3] = chiP0 + chiP1 + 2*(chi0 - chi1);
                CoeffsLifeOut[t*(aN+1)*4 + aIdx*4 + 0] = Coeffs[(aIdx*4) + 0];
                CoeffsLifeOut[t*(aN+1)*4 + aIdx*4 + 1] = Coeffs[(aIdx*4) + 1];
                CoeffsLifeOut[t*(aN+1)*4 + aIdx*4 + 2] = Coeffs[(aIdx*4) + 2];
                CoeffsLifeOut[t*(aN+1)*4 + aIdx*4 + 3] = Coeffs[(aIdx*4) + 3];
            }
            else { /* Calculate the exterior interval coefficients */
                Coeffs[0] = chiNow[0] - chiPNow[0]*muGrid[0];
                Coeffs[1] = chiPNow[0];
                Coeffs[2] = 0;
                Coeffs[3] = 0;
                CoeffsLifeOut[t*(aN+1)*4 + 0] = Coeffs[0];
                CoeffsLifeOut[t*(aN+1)*4 + 1] = Coeffs[1];
                CoeffsLifeOut[t*(aN+1)*4 + 2] = 0;
                CoeffsLifeOut[t*(aN+1)*4 + 3] = 0;
                Coeffs[aN*4 + 0] = chiNow[aN - 1] - chiPNow[aN - 1]*muGrid[aN - 1];
                Coeffs[aN*4 + 1] = chiPNow[aN - 1];
                Coeffs[aN*4 + 2] = 0;
                Coeffs[aN*4 + 3] = 0;
                CoeffsLifeOut[t*(aN+1)*4 + aN*4 + 0] = Coeffs[aN*4 + 0];
                CoeffsLifeOut[t*(aN+1)*4 + aN*4 + 1] = Coeffs[aN*4 + 1];
                CoeffsLifeOut[t*(aN+1)*4 + aN*4 + 2] = 0;
                CoeffsLifeOut[t*(aN+1)*4 + aN*4 + 3] = 0;
            }
        }
        barrier(CLK_GLOBAL_MEM_FENCE);
    } /* End of time loop */
}";


GPUkernelSolveLifeCycle = OpenCLFunctionLoad[LifeCycleSolveSrc
, "LifeCycleSolve_kernel"
, {
  {TypeWordB,1,"Output"}
, {TypeWordB,1,"Output"}
, {TypeWordB,1,"Input"}
, TypeWordB
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
, {TypeWordB,1,"Input"}
}
, {(NumOf\[Theta]ShockPoints+1)*NumOf\[CapitalPsi]ShockPoints*(\[GothicA]GridN + 1)}
(*, "Defines" -> {"USING_OPENCL_FUNCTION"->1,WhichGPUToUse(*,"mint"->"int","Real_t"->"float"*)}*)
(* This is part of the 'which GPU and which platform am I on' stuff that needs to be rethought *)
, Platform -> MyOpenCLPlatform
, Device -> MyOpenCLDevice
, "ShellOutputFunction"->Print];


(* MNW: This fills in the necessary information to pass to the GPU simulation kernel *)
GPUsimFunc[WealthVec_] := Block[{MoneyVec, ConsVec, AssetVec, NewWealthVec},

WealthMatrix = ConstantArray[0,{SimPeople*SimPeriods}];
MoneyMatrix = ConstantArray[0,{SimPeople*SimPeriods}];
ConsMatrix = ConstantArray[0,{SimPeople*SimPeriods}];
AssetMatrix = ConstantArray[0,{SimPeople*SimPeriods}];

{MoneyVec, ConsVec, AssetVec, NewWealthVec} = GPUkernelSim[
  MoneyMatrix, ConsMatrix, AssetMatrix, WealthMatrix
, WealthVec
, Flatten[\[Theta]SimMat]
, Flatten[\[CapitalPsi]SimMat]
, Flatten[Reverse[Take[\[Mu]GridsLife,-SimPeriods]]]
, Flatten[Reverse[Take[Map[Flatten[#] &, \[Chi]\[Mu]CoeffsLife],-SimPeriods]]]
, Reverse[Take[\[GothicA]LowerBoundLife,-SimPeriods]]
, Reverse[Take[\[FilledUpTriangle]\[GothicH]AccessibleLife,-SimPeriods]]
, Reverse[Take[\[Kappa]MinLife,-SimPeriods]]
, Take[\[ScriptCapitalG]Vect,SimPeriods]
];
(* 4 outputs, and 9 inputs *)

Return[{
  Partition[MoneyVec, SimPeople]
, Partition[ConsVec, SimPeople]
, Partition[AssetVec, SimPeople]
, Partition[NewWealthVec, SimPeople]

}];

ClearAll[WealthMatrix
, MoneyMatrix 
, ConsMatrix 
, AssetMatrix];

];


(* MNW: This is the source code for the GPU simulation routine. *)
NGridPoints = \[GothicA]GridN + 1 ;
(* Here NGridPoints=6 *)
(* It is not affected by whether 
we augment the grid at the far-left and the far-right end. *)

SimSrcMulti = "#ifdef USING_DOUBLE_PRECISIONQ
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif /* USING_DOUBLE_PRECISIONQ */

__kernel void SimMultiPeriod_kernel(
  __global " <> TypeWordA <> " * mtSimMatOut
, __global " <> TypeWordA <> " * ctSimMatOut
, __global " <> TypeWordA <> " * atSimMatOut
, __global " <> TypeWordA <> " * wtSimMatOut
, __global " <> TypeWordA <> " * WealthInit
, __global " <> TypeWordA <> " * TempShocks
, __global " <> TypeWordA <> " * PermShocks
, __global " <> TypeWordA <> " * muGrid
, __global " <> TypeWordA <> " * Coeffs
, __global " <> TypeWordA <> " * LowerBound
, __global " <> TypeWordA <> " * Deltah
, __global " <> TypeWordA <> " * PFMPCMin
, __global " <> TypeWordA <> " * Gamma) {

    int index = get_global_id(0);
    if (index >= " <> ToString[SimPeople] <>") return;
    " <> TypeWordA <> " Wealth = WealthInit[index];
    int t;
    " <> TypeWordA <> " Theta;
    " <> TypeWordA <> " Psi;
    " <> TypeWordA <> " Bound;
    " <> TypeWordA <> " Resources; 
    " <> TypeWordA <> " DeltaM;
    " <> TypeWordA <> " mu;
    " <> TypeWordA <> " muX;
    int SegmentStart;
    int ii;
    " <> TypeWordA <> " b0;
    " <> TypeWordA <> " b1;
    " <> TypeWordA <> " b2;
    " <> TypeWordA <> " b3;
    " <> TypeWordA <> " chi;
    " <> TypeWordA <> " Q;
    " <> TypeWordA <> " Cons;
    " <> TypeWordA <> " Assets;
    " <> TypeWordA <> " NextWealth;
    int Loc;
    int Temp;

    for (t = 0; t<" <> ToString[SimPeriods] <>"; t++) {
        Loc = t*" <> ToString[SimPeople] <>" + index;
        Theta = TempShocks[Loc];
        Psi = PermShocks[Loc];
        Bound = LowerBound[t];
        Resources = Wealth + Theta;
        DeltaM = Resources - Bound;
        mu = log(DeltaM);
        SegmentStart = t*(" <> ToString[NGridPoints] <>"+1)*4;
        Temp = 0;
        for (ii = 0; ii<" <> ToString[NGridPoints] <>"; ii++) {
            if (mu > muGrid[t*" <> ToString[NGridPoints] <>" + ii]) {
				SegmentStart += 4;
                Temp += 1;
            }
        };

        b0 = Coeffs[SegmentStart+0];
        b1 = Coeffs[SegmentStart+1];
        b2 = Coeffs[SegmentStart+2];
        b3 = Coeffs[SegmentStart+3];
        if ((Temp > 0) & (Temp < " <> ToString[NGridPoints] <>")) {
            muX = (mu - muGrid[t*" <> ToString[NGridPoints] <>" + Temp - 1])/(muGrid[t*" <> ToString[NGridPoints] <>" + Temp] - muGrid[t*" <> ToString[NGridPoints] <>" + Temp - 1]);
        }
        else {
            muX = mu;
        }

		chi = b0 + muX*(b1 + muX*(b2 + muX*(b3)));
        Q = 1/(1 + exp(chi));
        Cons = PFMPCMin[t]*(DeltaM + (1 - Q)*Deltah[t]);
		Cons = min(Cons, Resources);
        /* 
        This is to impose the 45 degree line constraint. 
        */

		Assets = Resources - Cons;
        NextWealth = Assets*(" <> ToString[\[DoubleStruckCapitalR]] <>"/(Gamma[t]*Psi));
        Wealth = NextWealth;

        mtSimMatOut[Loc]=Resources;
        ctSimMatOut[Loc]=Cons;
        atSimMatOut[Loc]=Assets;
        wtSimMatOut[Loc] = NextWealth;

        barrier(CLK_GLOBAL_MEM_FENCE);
    }
}";



(* MNW: This loads the GPU simulation kernel from the source code in SimSrcMulti. *)
(* CDC: For some machines with Nvidia GPUs, the kernel wrongly tries to compile assuming ATI; the line beginning Defines-> fixes this *)

GPUkernelSim = OpenCLFunctionLoad[SimSrcMulti
, "SimMultiPeriod_kernel"
, {
  {TypeWordB, 1, "Output"}
, {TypeWordB, 1, "Output"}
, {TypeWordB, 1, "Output"}
, {TypeWordB, 1, "Output"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
, {TypeWordB, 1, "Input"}
}
, {100} (* Workgroup size; experiments indicate 100 is a pretty good number for this, faster than much smaller or larger *)
(*, "Defines" -> {"USING_OPENCL_FUNCTION"->1,WhichGPUToUse(*"mint"->"int","Real_t"->"float"*)}*)
, Platform -> MyOpenCLPlatform
, Device -> MyOpenCLDevice
, "ShellOutputFunction"->Print];
(* 4 outputs, and 9 inputs *)
