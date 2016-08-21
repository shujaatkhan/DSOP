% ObjectFunction.m
% Function to be used in estimating param(s)

function F = ObjectFunction(x)

global PeriodsToSolve ...
    RFree n InitialWYRatio InitialWYRatioProb...
    VarInitialLogInc ProbOfAlive LevelAdjustingParameter...
    ThetaMat ThetaVec ThetaMatProb ThetaVecProb PermMat PermVec PermVecProb ...
    AlphaVec nP uP GList Betacorr NumOfPeriodsToSimulate NumOfPeople...
    ThetaList PermList stIndicator pi WealthCollege weight...
    RLife GammaLife lambdaMax KappaMin cLife DeltaGothicHLife ...
    GothicHMinLife GothicHExpLife GothicALowerBoundLife mLowerBoundLife ...
    chiIntData IntData Constrained

% Construct matrix of interpolation data
%ConstructcInterpFunc_EG
ConstructcInterpFunc_NoMoM
Simulate_NoMoM

F = WeightedSumDist(WealthCollege,stMedianListBy5Yrs,pi,weight);           % Evaluate function value


