function GothicV = GothicV(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Gothic V - The value function next period                                %
%                                                                          %
%   Inputs:                                                                %
%       a - beginning of period assets                                     %
%       Rho - Coefficient of Relative Risk Aversion                        %
%       Beta -  pure time discount factor                                  %
%       RFree - Interst Factor (1+r)                                       %
%       Gamma - Income growth                                              %
%       NumOfThetaShockPoints - Number of shocks in discrete approx.       %
%       ThetaVals - Value of transitory shocks   
%       PermVals - Value of permanent shocks
%   Outputs:                                                               %
%       GothicV - Value function next period                               %
%                                                                          %
%__________________________________________________________________________%
GothicV = zeros(size(a));

if PeriodsUntilT ~= 0  
for i = 1:length(ThetaVals)
    for j = 1:length(PermVals)
        mtp1 = (RLife(PeriodsUntilT)/(GammaLife(PeriodsUntilT)*PermVals(j)).*a + ThetaVals(i));
        GothicV = GothicV + (PermVals(j)*GammaLife(PeriodsUntilT))^(1-Rho)*PermVecProb(j)*ThetaVecProb(i)*...
            ScriptV(mtp1,Rho,BetaLife,RLife,GammaLife,ThetaVals,PermVals,PeriodsUntilT-1);
    %GothicV = GothicV + u(RFree * a + ThetaVals(j));
    end
end
GothicV = GothicV.*BetaLife(PeriodsUntilT);
end
% need to define ScriptV
% This is problematic if ThetaVals and PermVals are time-varying