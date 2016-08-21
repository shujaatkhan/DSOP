% Simulate.m
% This file runs simulation (slower version)

stList = zeros(NumOfPeriodsToSimulate,NumOfPeople);
mtList = zeros(NumOfPeriodsToSimulate,NumOfPeople);
ctList = zeros(NumOfPeriodsToSimulate,NumOfPeople);

% First period: 25
stList(1,:) = randsample(InitialWYRatio,NumOfPeople,true);                 % list of normalized s (savings at the beginning of age)      
mtList(1,:) = stList(1,:)+ThetaList(1,:);                                  % list of normalized m (cash on hand)

for i=1:NumOfPeople
    ctList(1,i) = interp1(IntData(:,1,end),IntData(:,2,end),mtList(1,i),'linear','extrap');  % consumption in period 25
%     if mtList(1,i)==0
%         ctList(1,i) = 0;
%     end
end

% Continue simulation
for t=2:NumOfPeriodsToSimulate                                             % age 26 onwards
    for j=1:NumOfPeople           
        stList(t,j) = (RLife(PeriodsToSolve-t+2)/GList(t-1)/PermList(t,j))*(mtList(t-1,j)-ctList(t-1,j));
        mtList(t,j) = stList(t,j) + ThetaList(t,j);
        ctList(t,j) = interp1(IntData(:,1,end-t+1),IntData(:,2,end-t+1),mtList(t,j),'linear','extrap');
    end
end
stMedianList      = (median(stList'))';                                    % list of median of savings level 
stTop25List = zeros(1,NumOfPeriodsToSimulate);
stBot25List = zeros(1,NumOfPeriodsToSimulate);

for i = 1:NumOfPeriodsToSimulate
stTop25List(i) = quantile(stList(i,:),0.75);
stBot25List(i) = quantile(stList(i,:),0.25);
end
% figure;
% plot(1:NumOfPeriodsToSimulate,stMedianList,'b',...
%     1:NumOfPeriodsToSimulate,stTop25List,'r--',...
%     1:NumOfPeriodsToSimulate,stBot25List,'k:')

for t=1:7
    stMedianListBy5Yrs(t) = mean(stMedianList((t-1)*5+2:t*5+1));  
    stTop25ListBy5Yrs(t) = mean(stTop25List((t-1)*5+2:t*5+1));
    stBot25ListBy5Yrs(t) = mean(stBot25List((t-1)*5+2:t*5+1));
end