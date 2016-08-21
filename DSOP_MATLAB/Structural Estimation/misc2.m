% setup_shocks.m

%% Construct the possible values of the shock to income
% Construct ThetaVec 
ThetaMat
ThetaMat(i,:)     = ThetaVec;

ThetaVecProb
ThetaVecRetireProb
ThetaMatProb = [repmat(ThetaVecProb,PeriodsToSolve-(90-65),1);repmat(ThetaVecRetireProb,90-65,1)];

PermMat = zeros(PeriodsToSolve,NumOfShockPoints);
PermVec     = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfShockPoints);
PermVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);

% clear unnecessary variables
clear ThetaVec ThetaVecProb PermVec

% SO left with:

ThetaMat
ThetaVecRetireProb
ThetaMatProb
PermMat
PermVecProb





