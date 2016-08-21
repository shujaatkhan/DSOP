function GothicCaa = GothicCaa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)

GothicCaa = (GothicCa(a+0.00001,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)-...
    GothicCa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT))./0.00001;