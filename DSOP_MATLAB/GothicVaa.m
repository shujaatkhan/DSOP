function GothicVaa = GothicVaa(a,Rho,BetaLife,RLife,GammaLife,mLowerBoundLife,DeltaGothicHLife,KappaMin,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT,chiIntData,Constrained)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% 2nd Derivative of Gothic V                                               %
%                                                                          %
%   Inputs:                                                                %
%       a - beginning of period assets                                     %
%       Rho - Coefficient of Relative Risk Aversion                        %
%       Beta -  pure time discount factor                                  %
%       RFree - Interst Factor (1+r)                                       %
%       NumOfThetaShockPoints - Number of shocks in discrete approx.       %
%       ThetaVals - Value of the shocks                                    %
%   Outputs:                                                               %
%       GothicVaa - 2nd Derivative of Gothic V                             %
%                                                                          %
%__________________________________________________________________________%

GothicVaa=zeros(size(a)); 

for i=1:length(ThetaVals)
    for j = 1:length(PermVals)
    mtp1 = (RLife(PeriodsUntilT+1)/(GammaLife(PeriodsUntilT+1)*PermVals(j)).*a + ThetaVals(i));
    GothicVaa = GothicVaa+((GammaLife(PeriodsUntilT+1)*PermVals(j))^(-Rho-1)).*...
        uPP(ScriptC(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained),Rho).*...
        ScriptKappa(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData)...
        *PermVecProb(j)*ThetaVecProb(i);
    end
end
GothicVaa = GothicVaa.*BetaLife(PeriodsUntilT+1)*RLife(PeriodsUntilT+1)^2;
end

% need to define function kappa