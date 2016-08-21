function GothicCa = GothicCa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)

GothicCa = GothicVaa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)./...
    uPP(GothicVa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)^(-1/Rho),Rho);