% setup_params.m
% Model parameters

% Set baseline values of model parameters
Constrained        = 0;                                                    % NOTE: Include this option
PeriodsToSolve     = 90-25;                                                % Number of periods to iterate
NumOfShockPoints   = 7;                                                    % Number of points in the discrete approximation to lognormal dist
pUnemp             = 0.5/100;                                              % Probability of unemployment (when unemployed inc level is zero)
pRetire            = 0+ind_mom*(pUnemp/10);                                % Probability of 'unemployment' after retirement (non-zero if MoM is used)
RFree              = 1.03;                                                 % Gross interest rate
AlphaMin           = 0.001;                                                % Minimum point in AlphaVec (glid of possible saving)
AlphaMax           = 20;                                                   % Maximum point in AlphaVec
AlphaHuge          = [];                                                   % Value of Alpha at which we connect to perf foresight function
n                  = 10;                                                   % Number of points in AlphaVec
%--------------------------------------------------------------------------
% Standard deviation of permanent and transitory shocks.
% There are no permanent shocks after retirement.
% There are no transitory shocks after retirement, unless MoM is used,
% which is when there is a small probability of a zero-income shock.
% NOTE: If MoM is used, the code produces an error if 0 st. dev. is used
% for any periods.
%--------------------------------------------------------------------------
if ind_mom && (pUnemp==0)
    error('If MoM is used, then pUnemp has to be non-zero.')
end
if ~ind_sdtv                                                               % Default
    SIGMAPerm = [0.1*ones(1,PeriodsToSolve-(90-65)) 0*ones(1,90-65)];      
    SIGMATran = [0.1*ones(1,PeriodsToSolve-(90-65)) 0*ones(1,90-65)];
else
    VarPerm = zeros(1,PeriodsToSolve-(90-65));
    VarTran = zeros(1,PeriodsToSolve-(90-65));
    %----------------------------------------------------------------------
    % Time-varying variance of shocks during working age, i.e. 26-65
    %----------------------------------------------------------------------
    for i = 1:length(VarPerm)
        if i<=20
            VarPerm(i) = 0.075-(0.075-0.03)/(PeriodsToSolve-(90-65)-20-1)*(i-1);
        else
            VarPerm(i) = 0.03+(0.08-0.03)/20*(i-20);
        end
        VarTran(i) = 0.14-(0.14-0.03)/(PeriodsToSolve-(90-65)-1)*(i-1);
    end
    SIGMAPerm = [sqrt(VarPerm) 0*ones(1,90-65)];
    SIGMATran = [sqrt(VarTran) 0*ones(1,90-65)];
end

InitialWYRatio     = [0.17, .5, .83];                                      % Initial wy ratio (from the program on the paper (p.13))
InitialWYRatioProb = [.33333, .33333, .33334];                             % Prob associated with initial wy ratio. (Using rand)
VarInitialLogInc   = .3329;                                                % Variance of initial log income

ProbOfAlive = [                                                            % Probability of being alive after retirement
 9.8438596e-01,9.8438596e-01,9.8438596e-01,9.8438596e-01,9.8438596e-01,... % (1st element is the prob of being alive until age 66)
 9.7567062e-01,9.7567062e-01,9.7567062e-01,9.7567062e-01,9.7567062e-01,...
 9.6207901e-01,9.6207901e-01,9.6207901e-01,9.6207901e-01,9.6207901e-01,...
 9.3721595e-01,9.3721595e-01,9.3721595e-01,9.3721595e-01,9.3721595e-01,...
 6.3095734e-01,6.3095734e-01,6.3095734e-01,6.3095734e-01,6.3095734e-01];
ProbOfAlive = [ones(1,PeriodsToSolve-length(ProbOfAlive)),ProbOfAlive];

pi = 0.50;                                                                 % Position of match (pith quantile is matched)
