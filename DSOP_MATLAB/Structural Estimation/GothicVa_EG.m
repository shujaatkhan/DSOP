function GothicVa = GothicVa_EG(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT,futureM, futureC,Constrained)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% This version does not employ the method of moments.                      %
% SAK: Update 08/09/2016                                                   %
%      This version does not require the mLowerBoundLife argument.         %
%                                                                          %
% Derivative of Gothic V - Derivative of the value function next period    %
%                                                                          %
%   Inputs:                                                                %
%       a - beginning of period assets                                     %
%       Rho - Coefficient of Relative Risk Aversion                        %
%       Beta -  pure time discount factor                                  %
%       RFree - Interst Factor (1+r)                                       %
%       Gamma - Income Growth                                              %
%       NumOfThetaShockPoints - Number of shocks in discrete approx.       %
%       ThetaVals - Value of the shocks                                    %
%   Outputs:                                                               %
%       GothicVa - Marginal Value function next period                     %
%                                                                          %
%__________________________________________________________________________%

GothicVa=zeros(size(a)); 

for i=1:length(ThetaVals)
    for j = 1:length(PermVals)
    
    mtp1 = (RLife(PeriodsUntilT+1)/(GammaLife(PeriodsUntilT+1)*PermVals(j)).*a + ThetaVals(i));
    
    %uP(ScriptC_EG(mtp1, mLowerBoundLife, futureM, futureC, Constrained),Rho)
    
    % Find the next-period consumption values at the mtp1 points:
    scriptc = interp1(futureM,futureC,mtp1,'linear','extrap');
    
    % Vectorized version of the constraint from ScriptC.m
    % TODO:
    % NOTE that there needs to be a warning here -- nothing *here* prevents us
    % from setting consumption negative for negative m-points. This should 
    % be handled by the selection of the m-space, but we need to be careful 
    % about this. This may be the source of the imaginary numbers seen when 
    % running constrained version of the Matlab code. 
    if Constrained == 1
        scriptc = min(scriptc, mtp1);  % NOTE: This is incorrect if mtp1 < 0, 
        % since we should never have negative consumption. Presumably we need to
        % have a line like the following; this was commented out in Script.m: 
        scriptc = max(scriptc,10^-6); % NOTE: commented out in ScriptC.m 
    end

    GothicVa = GothicVa+((GammaLife(PeriodsUntilT+1)*PermVals(j))^(-Rho))*...
                uP(scriptc,Rho)...
                *PermVecProb(j)*ThetaVecProb(i);
    end
end

GothicVa = GothicVa.*BetaLife(PeriodsUntilT+1)*RLife(PeriodsUntilT+1);
end


% need to define ScriptC, like Cnextp



