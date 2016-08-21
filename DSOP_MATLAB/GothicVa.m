function GothicVa = GothicVa(a,Rho,BetaLife,RLife,GammaLife,mLowerBoundLife,DeltaGothicHLife,KappaMin,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT,chiIntData,Constrained)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
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
    GothicVa = GothicVa+((GammaLife(PeriodsUntilT+1)*PermVals(j))^(-Rho))*...
                uP(ScriptC(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained),Rho)...
        *PermVecProb(j)*ThetaVecProb(i);
    end
end

GothicVa = GothicVa.*BetaLife(PeriodsUntilT+1)*RLife(PeriodsUntilT+1);
end


% need to define ScriptC, like Cnextp



