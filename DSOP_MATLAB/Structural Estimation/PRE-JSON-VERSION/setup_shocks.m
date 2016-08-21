% setup_shocks.m

%% Construct the possible values of the shock to income
% Construct ThetaVec
ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints);
if pUnemp>0
    ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints+1);
    % Add additional space for 0-valued shock if needed.
end

% Loop to fill in the ThetaMat matrix for temporary shocks. These will be the
% nodes (i.e. points of the discretized pmf) at which we will evaluate
% expectations. The associated probabilites (i.e. weights associated with each
% node) will be defined below.
%
% Note: the retirement and non-retirement ages all face the same style of
% uncertainty
for i = 1:PeriodsToSolve

 % Create a temporary holder of discrete pmf points; 1xn row vector:
 ThetaVec     = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfShockPoints);

 % If we assume unemployment, we need to add the "unemployment" value here, and
 % adjust all other values to preserve the expected value:
 if pUnemp>0
     if i> PeriodsToSolve-(90-65)
         ThetaVec = [0,ThetaVec/(1-pRetire)];  % During retirement years
     else
         ThetaVec = [0,ThetaVec/(1-pUnemp)];   % During regular years
     end
 end
ThetaMat(i,:)     = ThetaVec;   % Save the temporary Theta nodes.
end

% Create the probabilites associated with the above values. Note that these
% will be the same across the entire "employment" set of poeriods, and across
% the "retirement" periods (although will likely differ between those epochs).
ThetaVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);

 % If we include the possibility of "unemployment," i.e. a zero-income shock,
 % then we need to also be sure to adjust the probabilities associated with the
 % nodes above.
 if pUnemp>0 % If assume unemployment
    ThetaVecProb = [pUnemp,ThetaVecProb*(1-pUnemp)];
 end

% Do the same for retirement probabilites:
ThetaVecRetireProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);
 if pUnemp>0 % If assume unemployment
    ThetaVecRetireProb = [pRetire,ThetaVecRetireProb*(1-pRetire)];
 end

% Repeat the retirement and non-retirement probabilites the correct number of
% times, to match up with the nodes created for ThetaMat above.
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
