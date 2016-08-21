function scriptv = ScriptV(m,Rho,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained)


if PeriodsUntilT == 0;
    scriptv = u(m,Rho);
else
    vfromX =scriptVfromX(m,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained);
    scriptv = vfromX;
end
