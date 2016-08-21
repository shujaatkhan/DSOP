% setup_params.m
% Model parameters 

% Set baseline values of model parameters
Constrained            = 0;
PeriodsToSolve     = 90-25;                        % Number of periods to iterate
NumOfShockPoints   = 8;                            % Number of points in the discrete approximation to lognormal dist
VarPerm = zeros(1,PeriodsToSolve-(90-65));
VarTran = zeros(1,PeriodsToSolve-(90-65));
for i = 1:length(VarPerm)
    if i<=20
    VarPerm(i) = 0.075-(0.075-0.03)/(PeriodsToSolve-(90-65)-20-1)*(i-1);
    else
        VarPerm(i) = 0.03+(0.08-0.03)/20*(i-20);
    end
    VarTran(i) = 0.14-(0.14-0.03)/(PeriodsToSolve-(90-65)-1)*(i-1);
end
SIGMAPerm          = [sqrt(VarPerm) 0.001*ones(1,90-65)]; 
SIGMATran          = [sqrt(VarTran) 0.001*ones(1,90-65)];  
% SIGMAPerm          = [0.128*ones(1,PeriodsToSolve-(90-65)) 0.001*ones(1,90-65)];                        % Standard deviation of lognormal distribution of perm shocks
% SIGMATran          = [0.164*ones(1,PeriodsToSolve-(90-65)) 0.001*ones(1,90-65)];                        % Standard deviation of lognormal distribution of tran shocks
 % Note: these parameters are from Carroll (1992), a little more conservative than Carroll and Samwick (1997) 
pUnemp             = 0;%0.5/100;                      % Probability of unemployment (when unemployed inc level is zero) 
pRetire            = pUnemp/10;
RFree               = 1.03;                         % Gross interest rate
AlphaMin           = 0.00001;                      % Minimum point in AlphaVec (glid of possible saving)
AlphaMax           = 4;                            % Maximum point in AlphaVec
AlphaHuge          = 5;                         % Value of Alpha at which we connect to perf foresight function
n                  = 10;                           % Number of points in AlphaVec
% InitialWYRatio = exp(-2.794+1.784^2/2)*DiscreteApproxToMeanOneLogNormal(1.784,30);
% InitialWYRatioProb = 1/30*ones(1,30);
% InitialWYRatio = [-0.5, .5, 2];
% InitialWYRatioProb = [.33333, .33333, .333334];
InitialWYRatio     = [0.17, .5, .83]; %[0.3 0.98 2.17];             % Initial wy ratio (from the program on the paper (p.13))
InitialWYRatioProb = [.33333, .33333, .33334];    % Prob associated with initial wy ratio 
VarInitialLogInc   = .3329;                        % Variance of initial log income

% Probability of being alive after retirement 
% (1st element is the prob of being alive until age 66)
ProbOfAlive = [9.8438596e-01   9.8438596e-01   9.8438596e-01   9.8438596e-01   9.8438596e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01];
ProbOfAlive = [ones(1,PeriodsToSolve-length(ProbOfAlive)),ProbOfAlive];

% Position of match 
pi              = 0.50; % pith quantile is matched 

