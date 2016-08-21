% ConstructcInterpFunc.m 

% Construct matrix of interpolation data (JYao: using method of moderation)
% C = (0:n+1)'*ones(1,PeriodsToSolve+1); 
% M = (0:n+1)'*ones(1,PeriodsToSolve+1); 

% lambdaMax = 0;       % Marginal propensity to save
% kappaMin = 1;        % Marginal propensity to consume 
% cLife = 1;           
% DeltahLife = 0;      % Excess human wealth: difference between expected and minimum human wealth
% hMinLife = 0;        % Minimum human wealth at the end of period
% hExpLife = 0;        % Expected human wealth at the beginning of period excluding this period's income
% yExpPDV = 0;         % Expected human wealth including this period's income of 1, set to 0 before iteration 
% yMinPDV = 0;         % Minimum human wealth (including this period's minimum income),  set to 0 before iteration 
% aLowerBoundLife = 0;
% mLowerBoundLife = 0;

RLife = ones(1,PeriodsToSolve)*RFree;
BetaLife = zeros(1,PeriodsToSolve);
GammaLife = fliplr(GList);
lambdaMax = zeros(1,PeriodsToSolve+1);       % Marginal propensity to save
KappaMin = ones(1,PeriodsToSolve+1);        % Marginal propensity to consume 
cLife = ones(1,PeriodsToSolve+1);           
DeltaGothicHLife = zeros(1,PeriodsToSolve+1);      % Excess human wealth: difference between expected and minimum human wealth
GothicHMinLife = zeros(1,PeriodsToSolve+1);        % Minimum human wealth at the end of period
GothicHExpLife = zeros(1,PeriodsToSolve+1);        % Expected human wealth at the beginning of period excluding this period's income
yExpPDV = zeros(1,PeriodsToSolve);         % Expected human wealth including this period's income of 1, set to 0 before iteration 
yMinPDV = zeros(1,PeriodsToSolve);         % Minimum human wealth (including this period's minimum income),  set to 0 before iteration 
GothicALowerBoundLife = zeros(1,PeriodsToSolve+1);
mLowerBoundLife = zeros(1,PeriodsToSolve+1);
GothicAVec = AlphaVec;
chiIntData = zeros(length(GothicAVec)+1,3,PeriodsToSolve+1);   % Values of the last period do not matter
IntData = zeros(length(GothicAVec)+1,3,PeriodsToSolve+1);
IntData(:,:,1) = [(1:length(GothicAVec)+1)' (1:length(GothicAVec)+1)' ones(length(GothicAVec)+1,1)];  % [mVec cVec kappaVec]  
Rho = x(1);


for l=1:PeriodsToSolve
  P            = ProbOfAlive(length(ProbOfAlive)-l+1);  % Probability of being alive until next period 
  BetaLife(l)  = x(2)*Betacorr(length(Betacorr)-l+1)*P;
  
  ThetaVec     = ThetaMat(length(ThetaMatProb)-l+1,:);
  ThetaVecProb = ThetaMatProb(length(ThetaMatProb)-l+1,:);
  PermVec      = PermMat(length(ThetaMatProb)-l+1,:);
%   PermVecProb  = 1/length(PermVec)*ones(size(PermVec));
  
%   yExpPDV = hExpLife(end) + 1;
%   yMinPDV = hMinLife(end)+ min(ThetaVec);
  
  yExpPDV(l) = GothicHExpLife(l) + 1;
  yMinPDV(l) = GothicHMinLife(l)+ min(ThetaVec);
%   yMinPDV(l) = min(yMinPDV(l),bl);
  
  if Constrained == 1
%       GothicHMinLife(l+1) = min(ThetaVec)*min(PermVec)*GammaLife(l)/RLife(l);
     GothicHMinLife(l+1) = min(PermVec)*yMinPDV(l)*GammaLife(l)/RLife(l);
  else
     GothicHMinLife(l+1) = min(PermVec)*yMinPDV(l)*GammaLife(l)/RLife(l);
  end
  
  GothicHExpLife(l+1) = yExpPDV(l)*GammaLife(l)/RLife(l);                      
  DeltaGothicHLife(l+1) = GothicHExpLife(l+1) - GothicHMinLife(l+1);          
  lambdaMax(l+1)  = RLife(l)^(1/Rho-1)*BetaLife(l)^(1/Rho);
  KappaMin(l+1) = 1/(1+lambdaMax(l+1)/KappaMin(l));
%   cLife(l+1) = 1+cLife(l)*lambdaMax(l+1);  %vSum
  GothicALowerBoundLife(l+1) = -GothicHMinLife(l+1);
  mLowerBoundLife(l+1) = -GothicHMinLife(l+1);
  
  
  aLowerBoundt = GothicALowerBoundLife(l+1);
    if min(ThetaVec) > 0 && Constrained == 1;
%         mt = (mt -min(ThetaVals)*GammaLife(end)+0.001)./RLife(end);
        GothicAVect = GothicAVec + aLowerBoundt;
        GothicAVect = [0,GothicAVect];
        GothicAVect = sort(GothicAVect);
    else
        GothicAVect = sort([0.005 GothicAVec])+ aLowerBoundt;
        GothicAVect = sort(GothicAVect);
    end;  
    
    GothicVaVect = GothicVa(GothicAVect,Rho,BetaLife,RLife,GammaLife,mLowerBoundLife,DeltaGothicHLife,KappaMin,ThetaVec,ThetaVecProb,PermVec,PermVecProb,l-1,chiIntData,Constrained);
    GothicVaaVect = GothicVaa(GothicAVect,Rho,BetaLife,RLife,GammaLife,mLowerBoundLife,DeltaGothicHLife,KappaMin,ThetaVec,ThetaVecProb,PermVec,PermVecProb,l-1,chiIntData,Constrained);
    cVect = nP(GothicVaVect,Rho);
    caVect = GothicVaaVect./uPP(cVect,Rho);
    kappaVect= caVect./(1+caVect);
    mVect = GothicAVect+cVect;
    deltamVect = mVect-mLowerBoundLife(l+1);
    muVect = log(deltamVect);
    
    IntData(:,:,l+1) = [mVect' cVect' kappaVect'];
    
    
    
%  % Calculate ct from each grid point in AlphaVec
%   ChiVec = nP(P*GothVP(AlphaVec,x(1),l),x(1)); % Inverse Euler equation, P*GothVP(a,rho) is expected value given savings amount a
%   MuVec  = AlphaVec+ChiVec;
%   M(:,l+1)      = [0,MuVec]';                  % Matrix of interpolation data
%   C(:,l+1)      = [0,ChiVec]';                 % Matrix of interpolation data
%   
  % Calculations for the previous period
  
%   hMinLife = [hMinLife min(PermVec)*min(ThetaVec)/R];   
%   hExpLife = [hExpLife yExpPDV/R];                      
%   DeltahLife = [DeltahLife hExpLife(end) - hMinLife(end)];          
%   lambdaMax  = [lambdaMax (Rhat*x(2)*Betacorr(length(Betacorr)-l+1)*P)^(1/x(1))/Rhat];
%   kappaMin = [kappaMin 1/(1+lambdaMax(end)/kappaMin(end))];
%   cLife = [cLife 1+cLife*lambdaMax(end)];
%   aLowerBoundLife = [aLowerBoundLife -hMinLife(end)];
%   mLowerBoundLife = [mLowerBoundLife -hMinLife(end)];
  
%   hMinLife(l+1) = min(PermVec)*min(ThetaVec)/R;   
  
  
  
  kappat = KappaMin(l+1);
  deltah = DeltaGothicHLife(l+1);
  cVectRealst = cVect;
  cVectOptmst = kappat.*(deltamVect+deltah);
  cVectPestmst = kappat.*deltamVect;
  kappaVectRealst = kappaVect;
  kappaVectOptmst = kappat;
  kappaVectPestmst = kappat;
  koppaVals = (cVectOptmst-cVectRealst)./(cVectOptmst-cVectPestmst);
  koppamuVals = deltamVect.*(kappaVectOptmst-kappaVectRealst)./(kappat*deltah);
  chiVals = log((1./koppaVals)-1);
  chimuVals = koppamuVals./((koppaVals-1).*koppaVals);
  chiIntData(:,:,l+1) = [muVect' chiVals' chimuVals'];

%   hMinLife(l+1) = min(PermVec)*min(ThetaVec)/R;   
%   hExpLife(l+1) = yExpPDV/R;                      
%   DeltahLife(l+1) = hExpLife(l+1) - hMinLife(l+1);          
%   lambdaMax(l+1)  = (Rhat*x(2)*Betacorr(length(Betacorr)-l+1)*P)^(1/x(1))/Rhat;
%   kappaMin(l+1) = 1/(1+lambdaMax(l+1)/kappaMin(l));
%   cLife(l+1) = 1+cLife(l+1)*lambdaMax(l+1);
%   aLowerBoundLife(l+1) = -hMinLife(l+1);
%   mLowerBoundLife(l+1) = -hMinLife(l+1);
  
  
end;


