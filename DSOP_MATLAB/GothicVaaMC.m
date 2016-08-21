function GothicVaaMC = GothicVaaMC(a,Rho,BetaLife,RLife,GammaLife,NumOfThetaShockPoints,ThetaVals,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained,calRVals)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Derivative of Gothic V - Derivative of the value function next period    %
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
%       thorn -                                                            %
%       chiIntData - data for interpolation of chi                         %
%       Constrained - indicator of whether the individual is constrained   %
%       GothicCInterpData - data for interpolation of Gothic C             %
%       yExpPDV -                                                          %
%       calRVals - values for return on the risky asset                    %
%       VarSigma - Share of portfolio in the risky asset                   %
%   Outputs:                                                               %
%       GothicVaOpt - Marginal Value function next period at optimal share %
%                                                                          %
%__________________________________________________________________________%
for k = 1:length(a)
    VarSigmaOpt(k) = VarSigmaRaw(a(k),Rho,BetaLife,RLife,GammaLife,NumOfThetaShockPoints,ThetaVals,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained,calRVals);
end
GothicVaaMC = zeros(size(a));
for i = 1:length(calRVals)
    for j = 1:NumOfThetaShockPoints
        if PeriodsUntilT == 0
            GothicVaaMC = 0;
        elseif PeriodsUntilT == 1
            for k = 1:length(a)
                bbR = RLife(PeriodsUntilT+1)+ (calRVals(i)-RLife(PeriodsUntilT+1)).*VarSigmaOpt(k);
                mtp1 = a(k).*bbR./GammaLife(PeriodsUntilT+1)+ThetaVals(j);
                kappa = Kappamult(a(k),mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT-1,chiIntData,Constrained);
                GothicVaaMC(k) = GothicVaaMC(k)+BetaLife(PeriodsUntilT+1).*((GammaLife(PeriodsUntilT+1)).^(-Rho-1)).*(1/NumOfThetaShockPoints).*(1/length(calRVals)).*(bbR.^2).*uPP(mtp1,Rho).*kappa;
            end;
        else
            for k = 1:length(a)
                bbR = RLife(PeriodsUntilT+1)+ (calRVals(i)-RLife(PeriodsUntilT+1)).*VarSigmaOpt(k);
                mtp1 = a(k).*bbR./GammaLife(PeriodsUntilT+1)+ThetaVals(j);
                con =ScriptC(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT-1,chiIntData,Constrained);
                kappa = Kappamult(a(k),mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT-1,chiIntData,Constrained);
                GothicVaaMC(k) = GothicVaaMC(k) + BetaLife(PeriodsUntilT+1).*((GammaLife(PeriodsUntilT+1)).^(-Rho-1)).*(1/NumOfThetaShockPoints).*(1/length(calRVals)).*(bbR.^2).*uPP(con,Rho).*kappa;
            end;
        end
    end
end