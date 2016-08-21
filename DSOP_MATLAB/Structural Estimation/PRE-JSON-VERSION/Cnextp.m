% Cnextp.m
% Ct+1 function 

function c = Cnextp(m,t)
% Cnextp is constructed by interpolation to be the next-period consumption function Ct+1()
% t = 1 is the last period. 

global M C GothicHMinLife GothicHExpLife DeltaGothicHLife kappaMin mLowerBoundLife

mtp1 = M(:,t);  % data for the next-period consumption function
ctp1 = C(:,t);  % data for the next-period consumption function

GothicHMinLifetp1 = GothicHMinLife(t);
DeltaGothicHLifetp1 = DeltaGothicHLife(t);
GothicHExpLifetp1 = GothicHExpLife(t);
kappaMintp1 = kappaMin(t);

c = zeros(size(m));

muVec = log(mtp1(2:end)+GothicHMinLifetp1);
if DeltaGothicHLifetp1 == 0;
    chiVec = zeros(size(mtp1(2:end)));
else
    QVec = (kappaMintp1*(mtp1(2:end)+GothicHExpLifetp1)-ctp1(2:end))./(DeltaGothicHLifetp1*kappaMintp1);
%     QmuVec = (mtp1(2:end)+hExpLifetp1).*()
    chiVec = log(1./QVec-1);
end

ctRealst = scriptCfromChi(m,mLowerBoundLife,DeltaGothicHLife,kappaMin,t-1,chiIntData);
c = ctRealst;
% % extrapolate above maximal m
% iAbove = m >= mtp1(end);
% chislopeAbove  = (chiVec(end)-chiVec(end-1))/(muVec(end)-muVec(end-1));
% c(iAbove)   = kappaMintp1*(m(iAbove)+hExpLifetp1)-1./(1+exp(chiVec(end)+chislopeAbove*(log(m(iAbove)+hMinLifetp1)-muVec(end))))*DeltahLifetp1*kappaMintp1;
% 
% % extrapolate below minimal m
% iBelow = m <= mtp1(1);
% slopeBelow  = 1;
% c(iBelow)   = ctp1(1) + (m(iBelow)-mtp1(1))*slopeBelow;
% 
% % iBelow = m <= mtp1(1);
% % chislopeBelow  = (chiVec(2)-chiVec(1))/(muVec(2)-muVec(1));
% % c(iBelow)   = kappaMin*(m(iBelow)+hExpLife)-1./(1+exp(chislopeBelow*(log(m(iBelow)+hMinLife)-muVec(1))))*DeltahLife*kappaMin;
% 
% 
% % interpolate
% iInterp = ~(iAbove | iBelow);
% c(iInterp)  = interp1(mtp1,ctp1,m(iInterp));  

