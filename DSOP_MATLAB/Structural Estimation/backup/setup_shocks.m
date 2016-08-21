% setup_shocks.m

%% Construct the possible values of the shock to income
% Construct ThetaVec 
ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints);
if pUnemp>0
    ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints+1);
end

for i = 1:PeriodsToSolve
ThetaVec     = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfShockPoints);
 if pUnemp>0 % If assume unemployment
     if i> PeriodsToSolve-(90-65)
         ThetaVec = [0,ThetaVec/(1-pRetire)];
     else
         ThetaVec = [0,ThetaVec/(1-pUnemp)];
     end
 end
ThetaMat(i,:)     = ThetaVec;
end

ThetaVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);
 if pUnemp>0 % If assume unemployment
    ThetaVecProb = [pUnemp,ThetaVecProb*(1-pUnemp)];
 end
ThetaVecRetireProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);
 if pUnemp>0 % If assume unemployment
    ThetaVecRetireProb = [pRetire,ThetaVecRetireProb*(1-pRetire)];
 end

ThetaMatProb = [repmat(ThetaVecProb,PeriodsToSolve-(90-65),1);repmat(ThetaVecRetireProb,90-65,1)];


% Construct PermVec
PermMat = zeros(PeriodsToSolve,NumOfShockPoints);
for i = 1:PeriodsToSolve
PermVec     = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfShockPoints);
PermMat(i,:)= PermVec;
end

PermVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);

% clear unnecessary variables
clear ThetaVec ThetaVecProb PermVec