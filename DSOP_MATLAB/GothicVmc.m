function GothicVmc = GothicVmc(a,Rho,BetaLife,RLife,GammaLife,NumOfThetaShockPoints,ThetaVals,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,CapitalChiIntData,thorn,chiIntData,Constrained,GothicCInterpData,calRVals)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Gothic V - The value function next period                                %
%                                                                          %
%   Inputs:                                                                %
%       a - beginning of period assets                                     %
%       Rho - Coefficient of Relative Risk Aversion                        %
%       BetaLife -  pure time discount factor                              %
%       RLife - Interst Factor (1+r)                                       %
%       GammaLife - Income growth                                          %
%       NumOfThetaShockPoints - Number of shocks in discrete approx.       %
%       mLowerBoundLife - Lower bound for money balances                   %
%       DeltaGothicHLife - excess human wealth                             %
%       KappaMin - MPC                                                     %
%       PeriodsUntilT - period                                             %
%       CapitalChiIntData - data for interpolation of Capital Chi          %
%       thorn - Absolute patience factor                                   %
%       chiIntData - data for interpolation of chi                         %
%       Constrained - indicator of whether the individual is constrained   %
%       GothicCInterpData - data for interpolation of Gothic C             %
%       calRVals - values for return on the risky asset                    %
%       VarSigma - Share of portfolio in the risky asset                   %
%   Outputs:                                                               %
%       VarSigmaRaw - Marginal Value function next period                  %
%                                                                          %
%__________________________________________________________________________%
for k = 1:length(a)
    VarSigmaOpt(k) = VarSigmaRaw(a(k),Rho,BetaLife,RLife,GammaLife,NumOfThetaShockPoints,ThetaVals,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained,calRVals);
end            
GothicVmc = zeros(size(a));
for i = 1:length(calRVals)
    for j = 1:NumOfThetaShockPoints
        if PeriodsUntilT == 0
            GothicVmc = 0;
        elseif PeriodsUntilT == 1
            for k = 1:length(a)
                bbR = RLife(PeriodsUntilT+1)+ (calRVals(i)-RLife(PeriodsUntilT+1)).*VarSigmaOpt(k);
                mtp1 = a(k).*bbR./GammaLife(PeriodsUntilT+1)+ThetaVals(j);
                GothicVmc(k) = GothicVmc(k)+u(mtp1,Rho);
            end;
        else
            for k = 1:length(a)
                bbR = RLife(PeriodsUntilT+1)+ (calRVals(i)-RLife(PeriodsUntilT+1)).*VarSigmaOpt(k);
                mtp1 = a(k).*bbR./GammaLife(PeriodsUntilT+1)+ThetaVals(j);
                GothicVmc(k) = GothicVmc(k) + v(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT-1,CapitalChiIntData,thorn,Rho,chiIntData,Constrained,GothicCInterpData);
            end;
        end
    end
end
GothicVmc = GothicVmc.*BetaLife(PeriodsUntilT+1).*((GammaLife(PeriodsUntilT+1)).^(1-Rho)).*(1/NumOfThetaShockPoints).*(1/length(calRVals));