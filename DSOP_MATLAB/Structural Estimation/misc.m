% setup_params.m
Constrained
PeriodsToSolve
NumOfShockPoints
VarPerm = zeros(1,PeriodsToSolve-(90-65));
VarTran = zeros(1,PeriodsToSolve-(90-65));
SIGMAPerm          = [sqrt(VarPerm) 0.001*ones(1,90-65)]; 
SIGMATran          = [sqrt(VarTran) 0.001*ones(1,90-65)];  
pUnemp
pRetire
RFree
AlphaMin
AlphaMax
AlphaHug % Value of Alpha at which we connect to perf foresight function
n = 10;                           % Number of points in AlphaVec

InitialWYRatio     = [0.17, .5, .83]; %[0.3 0.98 2.17];             % Initial wy ratio (from the program on the paper (p.13))
InitialWYRatioProb = [.33333, .33333, .33334];    % Prob associated with initial wy ratio 
VarInitialLogInc   = .3329;                        % Variance of initial log income

ProbOfAlive = [ones(1,PeriodsToSolve-length(ProbOfAlive)),ProbOfAlive];
% Position of match 
pi              = 0.50; % pith quantile is matched 

% -------------------------------------------------------------setup_functions.m
uP = inline('c.^(-rho)','c','rho');   % CRRA marginal utility function
nP = inline('c.^(-1/rho)','c','rho'); % Inverse of the CRRA marginal 


% -----------------------------------------------------------------setup_grids.m
AlphaVec = exp(exp(exp(linspace(log(log(log(AlphaMin+1)+1)+1),log(log(log(AlphaMax+1)+1)+1),n))-1)-1)-1;
AlphaVec = [AlphaVec,AlphaHuge];




% ----------------------------------------------------------------setup_shocks.m

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



% ---------------------------------------------------EstimateRhoAndBeta setup...

% disp('sigma of Perm shocks, sigma of Tran shocks, prob of unemployment')
% disp([SIGMAPerm,SIGMATran,pUnemp])
setup_GList;                        % Load GList
setup_Betacorr;                     % Load corrected beta
NumOfPeriodsToSimulate = 90-55+1;   % Length of life in simulation (simulate until age 60)
NumOfPeople            = 10000;      % Number of people to simulate
NumOfPeopleBootstrap   = 10000;      % Number of people to simulate in bootstopping process
NumOfBootstrap         = 30;        % Number of times to iterate bootstapping process
setup_shocksLists                   % Set up shock lists



