function GothicCaaa = GothicCaaa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)

GothicCaaa = (GothicCaa(a+0.0000001,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)...
            -GothicCaa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT))./0.0000001;
        
end
