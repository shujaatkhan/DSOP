function scriptVfromX = scriptVfromX(mt,Rho,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained);

deltamt = mt-mLowerBoundLife(PeriodsUntilT+1);
deltaht = DeltaGothicHLife(PeriodsUntilT+1);
kappat  = KappaMin(PeriodsUntilT+1);
mut     = log(deltamt);
Xt      = XFuncLife(chiIntData,mut,PeriodsUntilT+1);
Ct      = cLife(PeriodsUntilT+1);
ctOptmst= kappat.*(deltamt+deltaht);
ctPestmst=kappat.* deltamt;
GammatOptmst = ctOptmst*(Ct)^(1/(1-Rho));
GammatPesmst = ctPestmst*(Ct)^(1/(1-Rho));
Koppat   = 1./(1+exp(Xt));
GammatRealst = GammatOptmst - Koppat*(GammatOptmst-GammatPesmst);    

scriptVfromX = u(GammatRealst,Rho);
