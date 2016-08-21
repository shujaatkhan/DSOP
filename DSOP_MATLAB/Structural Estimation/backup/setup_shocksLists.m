% setup_shocksLists.m

%% Construct income shock draw lists

%% Construct income shock lists
ThetaList = zeros(NumOfPeriodsToSimulate,NumOfPeople);
PermList = zeros(NumOfPeriodsToSimulate,NumOfPeople);

for i     = 1:NumOfPeriodsToSimulate
%     tTepm = randperm(NumOfPeople);
%     pTepm = randperm(NumOfPeople);
%     ThetaDraws = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfPeople);
%      if pUnemp>0 % If assume unemployment
%      ThetaDraws = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfPeople-pUnemp*NumOfPeople);
%      ThetaDraws = ThetaDraws./(1 - pUnemp);
%      ThetaDraws = [0*ones(1,pUnemp*NumOfPeople),ThetaDraws];
%      end
% 
% % Perm shock draw list
%      PermShockDraws = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfPeople);
%     for j=1:NumOfPeople
%         ThetaList(i,j) = ThetaDraws(tTepm(j));     % List of Theta (tran shock) 
%         PermList(i,j)  = PermShockDraws(pTepm(j)); % List of perm shock 
%     end
pTepm = NumOfShockPoints*rand(1,NumOfPeople);
    pTepmIndex = max(1,ceil(pTepm));
    
    if pUnemp == 0
        tTepm = NumOfShockPoints*rand(1,NumOfPeople);
        tTepmIndex = max(1,ceil(tTepm));
        ThetaDraws = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfShockPoints);
    elseif pUnemp>0 % If assume unemployment
        tTepm = NumOfShockPoints*rand(1,NumOfPeople)/(1 - pUnemp);
        tTepmIndex = max(1,ceil(tTepm));
        tTepmIndex((tTepmIndex>NumOfShockPoints))=NumOfShockPoints+1;
     ThetaDraws = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfShockPoints);
     ThetaDraws = ThetaDraws./(1 - pUnemp);
     ThetaDraws = [ThetaDraws 0];
    end

% Perm shock draw list
    PermShockDraws = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfShockPoints);
    for j=1:NumOfPeople
        ThetaList(i,j) = ThetaDraws(tTepmIndex(j));     % List of Theta (tran shock) 
        PermList(i,j)  = PermShockDraws(pTepmIndex(j)); % List of perm shock 
    end
end

%% Construct wtIndicator (list of indicators for initial wealth)
stIndicator = zeros(1,NumOfPeople);
for i=1:NumOfPeople
    r = rand;
    stIndicator(i) = length(InitialWYRatioProb)+1-sum((r<=cumsum(InitialWYRatioProb)));
%     if r < InitialWYRatioProb(1)
%         stIndicator(i) = 1;
%     elseif r < InitialWYRatioProb(1)+InitialWYRatioProb(2)
%         stIndicator(i) = 2;
%     else
%         stIndicator(i) = 3;
%     end
end