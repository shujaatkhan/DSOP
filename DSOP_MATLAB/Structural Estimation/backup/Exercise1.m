
clear all; 
% close all;
TimeS = cputime;
disp('=======================')

% Setup 
setup_everything
% disp('sigma of Perm shocks, sigma of Tran shocks, prob of unemployment')
% disp([SIGMAPerm,SIGMATran,pUnemp])
setup_GList;                        % Load GList
setup_Betacorr;                     % Load corrected beta
NumOfPeriodsToSimulate = 90-55+1;   % Length of life in simulation (simulate until age 60)
NumOfPeople            = 20000;     % Number of people to simulate
NumOfPeopleBootstrap   = 1000;      % Number of people to simulate in bootstopping process
NumOfBootstrap         = 30;        % Number of times to iterate bootstapping process
setup_shocksLists                   % Set up shock lists


% Estimation
Data_SCF_wealth                          % Load data
WealthCollege = WealthPopulationCollege; % Default WealthCollege is population itself


%x = [3.8423 0.9891];  % zero natural borrowing constraint (unemployment prob > 0)  
x = [2.6579    0.9341    0.9008];%[4.0255 0.9816];   % nonzero natural borrowing constraint (unemployment prob = 0)
% Construct matrix of interpolation data   
ConstructcInterpFunc              
% Simulate 
Simulate             


%% How model does compared to the data
age = 30:5:60;
figure;
plot(age,stMedianListBy5Yrs,age,MeanWealthSCF,'b-.',age,MedianWealthSCF,'k',...
age,Top25WealthSCF,'g',age, stTop25ListBy5Yrs,'g--',...
age,Bot25WealthSCF,'r',age, stBot25ListBy5Yrs,'r--');
legend('Simulated Median', 'Data Mean', 'Data Median',...
    'Data Top25','Simulated Top25',...
    'Data Bot25','Simulated Bot25','Location','NorthWest')
ylabel('W/Y')
xlabel('Age')



% scatter(WealthCollege(:,2),WealthCollege(:,1),'r');
% hold off;

%% Total Wealth to Income Ratio of the model economy
SimulateExercise1

TotalWYRatio = sum(sum(stList.*ptList))/sum(sum(ytList));

%% Consumption insurance over the life cycle
CtList = ctList.*ptList;
% Permanent shocks
phi_perm = zeros(1,NumOfPeriodsToSimulate-1);
for i = 1:NumOfPeriodsToSimulate-1
    temp = cov(log(CtList(i+1,:))-log(CtList(i,:)),log(PermList(i+1,:)))/var(log(PermList(i+1,:)));
    phi_perm(i) = 1-temp(1,2);
end

% Transitory shocks
phi_tran = zeros(1,NumOfPeriodsToSimulate-1);
ThetaListNoZero = max(ThetaList,0.001);
for i = 1:NumOfPeriodsToSimulate-1
    index = (ThetaListNoZero(i+1,:)~=0);
    temp = cov(log(CtList(i+1,index))-log(CtList(i,index)),log(ThetaListNoZero(i+1,index)))/var(log(ThetaListNoZero(i+1,index)));
    phi_tran(i) = 1-temp(1,2);
end

age = 26:1:60;
figure;
plot(age,phi_perm,age,phi_tran,'r')
title('Insurance Coefficient')
legend('Permanent','Transitory')

% Average insurance coefficient 
fprintf('---------------------------Average insurance coefficient----------------------\n')
fprintf('Permanent      Transitory\n')
fprintf('%1.4f          %1.4f\n', mean(phi_perm), mean(phi_tran))

    
%% BPP method estimation
phi_perm_bpp = zeros(1,NumOfPeriodsToSimulate-3);
for i = 3:NumOfPeriodsToSimulate-1
    index = (ytList(i,:)~=0)&(ytList(i-1,:)~=0)&(ytList(i-2,:)~=0)&(ytList(i+1,:)~=0);
    temp1 = cov(log(ytList(i,index))-log(ytList(i-1,index)),log(ytList(i+1,index))-log(ytList(i-2,index)));
    var_perm = temp1(1,2);
    temp2 = cov(log(CtList(i,index))-log(CtList(i-1,index)),log(ytList(i+1,index))-log(ytList(i-2,index)));
    cov_Cperm = temp2(1,2);
    phi_perm_bpp(i-2) = 1-cov_Cperm/var_perm;
end


phi_tran_bpp = zeros(1,NumOfPeriodsToSimulate-3);
for i = 3:NumOfPeriodsToSimulate-1
    index = (ytList(i,:)~=0)&(ytList(i-1,:)~=0)&(ytList(i+1,:)~=0);
    temp1 = cov(log(ytList(i,index))-log(ytList(i-1,index)),log(ytList(i+1,index))-log(ytList(i,index)));
    var_tran = temp1(1,2);
    temp2 = cov(log(CtList(i,index))-log(CtList(i-1,index)),log(ytList(i+1,index))-log(ytList(i,index)));
    cov_Ctran = temp2(1,2);
    phi_tran_bpp(i-2) = 1-cov_Ctran/var_tran;
end
age = 27:1:59;
figure;
plot(age,phi_perm(2:end-1),age,phi_tran(2:end-1),'r',age,phi_perm_bpp,'b--',age,phi_tran_bpp,'r--');
title('Insurance Coefficient')
legend('Permanent','Transitory','Perm BPP','Tran BPP','Location','NorthWest')

figure;
plot(age,phi_perm_bpp-phi_perm(2:end-1),age,phi_tran_bpp-phi_tran(2:end-1),'r');
title('BPP error')
legend('Permanent','Transitory')


% Average insurance coefficient 
fprintf('---------------------------Average insurance coefficient----------------------\n')
fprintf('               Permanent      Transitory \n')
fprintf('Model True:    %1.4f          %1.4f\n', mean(phi_perm), mean(phi_tran))
fprintf('Model BPP:     %1.4f          %1.4f\n', mean(phi_perm_bpp), mean(phi_tran_bpp))

%% Log linearized Euler equation
mpc_perm = zeros(1,NumOfPeriodsToSimulate-1);
mpc_tran = zeros(1,NumOfPeriodsToSimulate-1);
for i = 1:NumOfPeriodsToSimulate-1
    deltaC = log(CtList(i+1,:))-log(CtList(i,:));
    permShock = log(PermList(i+1,:));
    tranShock = log(ThetaList(i+1,:));
    coeff = lscov([permShock' tranShock'],deltaC');
    mpc_perm(i) = coeff(1);
    mpc_tran(i) = coeff(2);
end
phi_perm_reg = 1- mpc_perm;
phi_tran_reg = 1- mpc_tran;

age = 27:1:59;
BPP = figure;
subplot(2,1,1);
% plot(age,phi_perm(2:end-1),age,phi_tran(2:end-1),'r',age,phi_perm_bpp,'b--',age,phi_tran_bpp,'r--',age,phi_perm_reg(2:end-1),'b-*',age,phi_tran_reg(2:end-1),'r-*');
% title('Insurance Coefficient')
% legend('Permanent','Transitory','Perm BPP','Tran BPP','Perm Reg','Tran Reg','Location','NorthWest')

plot(age,phi_perm(2:end-1),age,phi_tran(2:end-1),'r',age,phi_perm_reg(2:end-1),'b-*',age,phi_tran_reg(2:end-1),'r-*');
title('Insurance Coefficient')
legend('Permanent','Transitory','Perm BPP','Tran BPP','Location','NorthWest')
xlabel('Age')

age = 26:1:60;
subplot(2,1,2)
plot(age,100*(phi_perm_reg-phi_perm)./phi_perm,age,100*(phi_tran_reg-phi_tran)./phi_tran,'r');
title('Percentage Deviation from True Insurance Coefficients')
legend('Permanent','Transitory')
xlabel('Age')

FigureExFilePath = fullfile(pwd,'../Figures/BPP.pdf');
print ('-dpdf',BPP,FigureExFilePath);

FigureExFilePath = fullfile(pwd,'../Figures/BPP.eps');
print ('-depsc',BPP,FigureExFilePath);

%-----------------Partial Insurance Parameter------------------------%
age = 27:1:59;
Partial = figure;
subplot(2,1,1);
% plot(age,phi_perm(2:end-1),age,phi_tran(2:end-1),'r',age,phi_perm_bpp,'b--',age,phi_tran_bpp,'r--',age,phi_perm_reg(2:end-1),'b-*',age,phi_tran_reg(2:end-1),'r-*');
% title('Insurance Coefficient')
% legend('Permanent','Transitory','Perm BPP','Tran BPP','Perm Reg','Tran Reg','Location','NorthWest')

plot(age,1-phi_perm(2:end-1),age,1-phi_tran(2:end-1),'r',age,1-phi_perm_reg(2:end-1),'b-*',age,1-phi_tran_reg(2:end-1),'r-*');
title('Partial Insurance Parameters')
xlabel('Age')
legend('Perm True','Tran True','Perm BPP','Tran BPP')

age = 26:1:60;
subplot(2,1,2)
plot(age,100*((1-phi_perm_reg)-(1-phi_perm))./(1-phi_perm),age,100*((1-phi_tran_reg)-(1-phi_tran))./(1-phi_tran),'r');
title('Percentage Deviation from True Partial Insurance Parameters')
xlabel('Age')
legend('Permanent','Transitory','Location','NorthWest')

FigureExFilePath = fullfile(pwd,'../Figures/Partial.pdf');
print ('-dpdf',Partial,FigureExFilePath);

FigureExFilePath = fullfile(pwd,'../Figures/Partial.eps');
print ('-depsc',Partial,FigureExFilePath);
%------------------End of Partial Insurance Parameter---------------------%
% Proposed estimator
mpc_perm_prop = zeros(1,NumOfPeriodsToSimulate-1);
mpc_tran_prop = zeros(1,NumOfPeriodsToSimulate-1);
for i = 1:NumOfPeriodsToSimulate-1
    deltaC = (log(CtList(i+1,:))-log(CtList(i,:)))./stList(i,:);
    permShock = log(PermList(i+1,:))./stList(i,:);
    tranShock = log(ThetaList(i+1,:))./stList(i,:);
    coeff = lscov([permShock' tranShock'],deltaC');
    mpc_perm_prop(i) = coeff(1);
    mpc_tran_prop(i) = coeff(2);
end
phi_perm_reg_prop = 1- mpc_perm_prop;
phi_tran_reg_prop = 1- mpc_tran_prop;

age = 27:1:59;
figure;
plot(age,phi_perm(2:end-1),age,phi_tran(2:end-1),'r',age,phi_perm_bpp,'b--',age,phi_tran_bpp,'r--',age,phi_perm_reg_prop(2:end-1),'b-*',age,phi_tran_reg_prop(2:end-1),'r-*');
title('Insurance Coefficient')
legend('Permanent','Transitory','Perm BPP','Tran BPP','Perm Reg Prop','Tran Reg Prop','Location','NorthWest')

age = 26:1:60;
figure;
plot(age,(phi_perm_reg_prop-phi_perm)./phi_perm,age,(phi_tran_reg_prop-phi_tran)./phi_tran,'r');
title('BPP regression error')
legend('Permanent','Transitory')


%% Wealth to income ratio distribution
age = 30;
figure;
hist(stList(age-25,:),1000);

[f, xi] = ksdensity(mtList(age-25,:));
figure()
plot(xi,f);


percentile = 0:2:100;
mt_percentile = prctile(mtList(age-25,:),percentile);  % use market wealth to differentiate
phi_perm = zeros(1,length(percentile)-1);
phi_perm_bpp = zeros(1,length(percentile)-1);
phi_perm_bpp_reg = zeros(1,length(percentile)-1);
for i = 2:length(percentile)
    index = (mtList(age-25,:)<mt_percentile(i))&(mtList(age-25,:)>=mt_percentile(i-1));
    temp = cov(log(CtList(age-25,index))-log(CtList(age-25-1,index)),log(PermList(age-25,index)))/var(log(PermList(age-25,index)));
    phi_perm(i-1) = 1-temp(1,2);
    temp1 = cov(log(ytList(age-25,index))-log(ytList(age-25-1,index)),log(ytList(age-25+1,index))-log(ytList(age-25-2,index)));
    var_perm = temp1(1,2);
    temp2 = cov(log(CtList(age-25,index))-log(CtList(age-25-1,index)),log(ytList(age-25+1,index))-log(ytList(age-25-2,index)));
    cov_Cperm = temp2(1,2);
    phi_perm_bpp(i-1) = 1-cov_Cperm/var_perm;
    deltaC = log(CtList(age-25,index))-log(CtList(age-25-1,index));
    permShock = log(PermList(age-25,index));
    tranShock = log(ThetaList(age-25,index));
    coeff = lscov([permShock' tranShock'],deltaC');
    phi_perm_bpp_reg(i-1) = 1-coeff(1);
end
figure;
plot(mt_percentile(2:end),phi_perm,mt_percentile(2:end),phi_perm_bpp,'r',mt_percentile(2:end),phi_perm_bpp_reg,'g');
legend('Model True','Model BPP')

age = 30;
percentile = 0:2:100;
at_percentile = prctile(mtList(age-25-1,:)-ctList(age-25-1,:),percentile);  % use beginning-period-of-asset to differential
phi_perm = zeros(1,length(percentile)-1);
phi_perm_bpp = zeros(1,length(percentile)-1);
for i = 2:length(percentile)
    index = (mtList(age-25-1,:)-ctList(age-25-1,:)<at_percentile(i))&(mtList(age-25-1,:)-ctList(age-25-1,:)>=at_percentile(i-1));
    temp = cov(log(CtList(age-25,index))-log(CtList(age-25-1,index)),log(PermList(age-25,index)))/var(log(PermList(age-25,index)));
    phi_perm(i-1) = 1-temp(1,2);
    temp1 = cov(log(ytList(age-25,index))-log(ytList(age-25-1,index)),log(ytList(age-25+1,index))-log(ytList(age-25-2,index)));
    var_perm = temp1(1,2);
    temp2 = cov(log(CtList(age-25,index))-log(CtList(age-25-1,index)),log(ytList(age-25+1,index))-log(ytList(age-25-2,index)));
    cov_Cperm = temp2(1,2);
    phi_perm_bpp(i-1) = 1-cov_Cperm/var_perm;
    deltaC = log(CtList(age-25,index))-log(CtList(age-25-1,index));
    permShock = log(PermList(age-25,index));
    tranShock = log(ThetaList(age-25,index));
    coeff = lscov([permShock' tranShock'],deltaC');
    phi_perm_bpp_reg(i-1) = 1-coeff(1);
end
BPPcrosssection = figure;
% plot(at_percentile(2:end),phi_perm,at_percentile(2:end),phi_perm_bpp_reg,'g',at_percentile(2:end),phi_perm_bpp,'r');
% legend('Model True','Reg BPP','Model BPP')

plot(at_percentile(2:end),phi_perm,at_percentile(2:end),phi_perm_bpp_reg,'g');
legend('Model True','BPP reg')
ylabel('Insurance Coefficient')
xlabel('Wealth/Income')
h = sprintf('Insurance Coefficient of Permanent Shocks at Age %d',age);
title(h);
FigureExFilePath = fullfile(pwd,'../Figures/BPPcrosssection.pdf');
print ('-dpdf',BPPcrosssection,FigureExFilePath);
FigureExFilePath = fullfile(pwd,'../Figures/BPPcrosssection.eps');
print ('-depsc',BPPcrosssection,FigureExFilePath);

%% Carroll 2009 JME
at = at_percentile;
MPCP = zeros(size(at));
InsuCoeff = zeros(size(at));  % evaluated at mean
InsuCoeffTrue = zeros(size(at)); 
age = 30-25;
for i = 1:length(at)
    MPCP(i) = ScriptC(RLife(PeriodsToSolve-age+2)*at(i)+0.5,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData,Constrained)...
        -ScriptKappa(RLife(PeriodsToSolve-age+2)*at(i)+0.5,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData)*RLife(PeriodsToSolve-age+2)*at(i);
    InsuCoeff(i) = ScriptKappa(RLife(PeriodsToSolve-age+1)*at(i)+1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData)*RLife(PeriodsToSolve-age+2)*at(i)...
        /ScriptC(RLife(PeriodsToSolve-age+1)*at(i)+1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData,Constrained);
    for j = 1:length(PermVec)
        for k = 1:length(ThetaVec)
            InsuCoeffTrue(i) = InsuCoeffTrue(i)+...
                ScriptKappa(RLife(PeriodsToSolve-age+1)*at(i)/PermVec(j)+ThetaVec(k),mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData)*RLife(PeriodsToSolve-age+2)*at(i)/PermVec(j)...
        /ScriptC(RLife(PeriodsToSolve-age+1)*at(i)/PermVec(j)+ThetaVec(k),mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsToSolve-age+1,chiIntData,Constrained)*1/length(PermVec)*1/length(ThetaVec);
        end
    end    
end
Carrollcrosssection = figure;
plot(at(2:end),phi_perm,at(2:end),InsuCoeff(2:end),'r',at_percentile(2:end),phi_perm_bpp_reg,'g')
title(sprintf('Insurance Coefficient of Permanent Shocks at Age %d', age+25));
ylabel('Insurance Coefficient')
xlabel('Wealth/Income')
legend('Model True', 'New Approx','BPP Approx')
FigureExFilePath = fullfile(pwd,'../Figures/Carrollcrosssection.pdf');
print ('-dpdf',Carrollcrosssection,FigureExFilePath);
FigureExFilePath = fullfile(pwd,'../Figures/Carrollcrosssection.eps');
print ('-depsc',Carrollcrosssection,FigureExFilePath);

figure;
plot(at(2:end),phi_perm,at(2:end),InsuCoeff(2:end),'g',at(2:end),InsuCoeffTrue(2:end),'r')

%% Age inequality profile of consumption
Age = 26:1:60;
VarLogConsumption = zeros(size(Age));
VarLogIncome = zeros(size(Age));
for i = 1:length(Age)
    VarLogConsumption(i) = var(log(CtList(i,:)));
    VarLogIncome(i) = var(log(ytList(i,:)));
end
figure;
plot(Age,VarLogConsumption,Age,VarLogIncome,'r')
title('Age-Inequality Profile of Consumption and Income')
ylabel('Cross-sectional Variance')
xlabel('Age')
legend('Consumption','Income')

%% Age inequality profile of wealth
Age = 26:1:60;
VarWealth = zeros(size(Age));
for i = 1:length(Age)
    VarWealth(i) = var(stList(i,:));
end
figure;
plot(Age,VarWealth)
title('Age-Inequality Profile of Wealth')
ylabel('Cross-sectional Variance of Wealth')
xlabel('Age')


StList = stList.*ptList;
VarLogWealth = zeros(size(Age));
InterQuatileWealth = zeros(size(Age));
InterQuatileWealthGroup = zeros(1,7);
InterQuatileWealthData = zeros(1,7);
for i = 1:length(Age)
    VarLogWealth(i) = var(log(StList(i,:)));
end
for i = 1:length(Age)+1
InterQuatileWealth(i) = iqr(stList(i,:));
end
for t = 1:7
    InterQuatileWealthGroup(t) = mean(InterQuatileWealth((t-1)*5+2:t*5+1));
    InterQuatileWealthData(t) = Top25WealthSCF(1,t)-Bot25WealthSCF(1,t);
end
figure;
plot(Age,VarLogWealth)
title('Age-Inequality Profile of Wealth')
ylabel('Cross-sectional Variance of Log Wealth')
xlabel('Age')


Match = figure;
age = 30:5:60;
subplot(2,1,1);
plot(age,stMedianListBy5Yrs,age,MeanWealthSCF,'b-.',age,MedianWealthSCF,'k',...
age,Top25WealthSCF,'g',age, stTop25ListBy5Yrs,'g--',...
age,Bot25WealthSCF,'r',age, stBot25ListBy5Yrs,'r--');
legend('Simulated Median', 'Data Mean', 'Data Median',...
    'Data Top25','Simulated Top25',...
    'Data Bot25','Simulated Bot25','Location','NorthWest')
ylabel('W/Y')
xlabel('Age')

Group = 1:1:7;
subplot(2,1,2);
plot(Group,InterQuatileWealthGroup,Group,InterQuatileWealthData,'r')
title('Interquatile Range of Wealth Distribution by Age')
ylabel('Interquatile Range')
xlabel('Age')
legend('Simulated','Data','Location','NorthWest')

FigureExFilePath = fullfile(pwd,'../Figures/Match.pdf');
print ('-dpdf',Match,FigureExFilePath);

FigureExFilePath = fullfile(pwd,'../Figures/Match.eps');
print ('-depsc',Match,FigureExFilePath);

VarLogWealthData(1) = (Wealth26_30).^2*WealthDens26_30'-((Wealth26_30)*WealthDens26_30')^2;
VarLogWealthData(2) = (Wealth31_35).^2*WealthDens31_35'-((Wealth31_35)*WealthDens31_35')^2;
VarLogWealthData(3) = (Wealth36_40).^2*WealthDens36_40'-((Wealth36_40)*WealthDens36_40')^2;
VarLogWealthData(4) = (Wealth41_45).^2*WealthDens41_45'-((Wealth41_45)*WealthDens41_45')^2;
VarLogWealthData(5) = (Wealth46_50).^2*WealthDens46_50'-((Wealth46_50)*WealthDens46_50')^2;
VarLogWealthData(6) = (Wealth51_55).^2*WealthDens51_55'-((Wealth51_55)*WealthDens51_55')^2;
VarLogWealthData(7) = (Wealth56_60).^2*WealthDens56_60'-((Wealth56_60)*WealthDens56_60')^2;
figure;
plot(Group,VarLogWealthData)


% VarLogWealthData(1) = (log(max(0.01,Wealth26_30))).^2*WealthDens26_30'-(log(max(0.01,Wealth26_30))*WealthDens26_30')^2;
% VarLogWealthData(2) = (log(max(0.01,Wealth31_35))).^2*WealthDens31_35'-(log(max(0.01,Wealth31_35))*WealthDens31_35')^2;
% VarLogWealthData(3) = (log(max(0.01,Wealth36_40))).^2*WealthDens36_40'-(log(max(0.01,Wealth36_40))*WealthDens36_40')^2;
% VarLogWealthData(4) = (log(max(0.01,Wealth41_45))).^2*WealthDens41_45'-(log(max(0.01,Wealth41_45))*WealthDens41_45')^2;
% VarLogWealthData(5) = (log(max(0.01,Wealth46_50))).^2*WealthDens46_50'-(log(max(0.01,Wealth46_50))*WealthDens46_50')^2;
% VarLogWealthData(6) = (log(max(0.01,Wealth51_55))).^2*WealthDens51_55'-(log(max(0.01,Wealth51_55))*WealthDens51_55')^2;
% VarLogWealthData(7) = (log(max(0.01,Wealth56_60))).^2*WealthDens56_60'-(log(max(0.01,Wealth56_60))*WealthDens56_60')^2;
% figure;
% plot(Group,VarLogWealthData)
% MtList = mtList.*ptList;
% VarLogWealth = zeros(size(Age));
% for i = 1:length(Age)
%     VarLogWealth(i) = var(log(MtList(i,:)));
% end
% figure;
% plot(Age,VarLogWealth)
% title('Age-Inequality Profile of Wealth')
% ylabel('Cross-sectional Variance of Log Wealth')
% xlabel('Age')

