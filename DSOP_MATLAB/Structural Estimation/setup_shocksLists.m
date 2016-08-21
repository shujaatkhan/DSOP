% setup_shocksLists.m
% Construct income shock draw lists

ThetaList = zeros(NumOfPeriodsToSimulate,NumOfPeople);
PermList = zeros(NumOfPeriodsToSimulate,NumOfPeople);

%--------------------------------------------------------------------------
% SAK 08/11/16
% NOTE: The shocks are for periods s=26 to 90 (inclusive); however, in the
% original code, we were also simulating consumption for s=25, using shocks
% for s=26. This resulted in each period's simulation using the next
% period's (incorrect) shock.
% I have modified the code, so that now the income shock draws are for s=25
% to s=60 (inclusive) and in s=25, the shock draws are equal to 1 for all
% the people being simulated. Starting s=26, shocks are drawn from the
% LogNormal distribution.
%
% In the original code, in each period, the shock draws for the people
% being simulated were drawn from just the shock points. I have changed
% this, so that there are n=NumOfPeople different shocks distributed evenly
% across the distribution. This is how it is in the notes.
%
% NOTE: On the website version of the Matlab code, since the st. dev. of
% income shocks was not time-varying, shocks were drawn from the
% distribution only once and they were allotted randomly to different
% people in each period. This made the program fast.
% In Matt's version of the code, draws were only made from the shock
% points, which also made the program very fast.
% In this case, when we allow for age-varying st. dev. of income shocks and
% we also draw n=NumOfPeople shocks for each age, the program become much
% slower.
% For NumOfPeople=1000:
% running setup_shocksLists.m w/ age-varying st. dev. takes about 0.74mins.
% running setup_shocksLists.m with constant st. dev. takes about  0.03mins.
%--------------------------------------------------------------------------
ThetaList = zeros(NumOfPeriodsToSimulate,NumOfPeople);
PermList = zeros(NumOfPeriodsToSimulate,NumOfPeople);

if pUnemp>0
    NumOfUnemp = ceil(pUnemp*NumOfPeople);                                 % Number of unemployed people before retirement
    NumOfUnempRet = ceil(pRetire*NumOfPeople);                             % Number of zero-income shock people after retirement
else
    NumOfUnemp = 0;
    NumOfUnempRet = 0;
end

if ~ind_sdtv                                                               % If st. dev. of income shocks is constant, then simulation becomes fast
    if NumOfPeriodsToSimulate>(65-25+1)
        error('Simulating retirement years as well. Write the code for retirement year draws.')
    end
    ThetaDraws_tmp = DiscreteApproxToMeanOneLogNormal(SIGMATran(1),NumOfPeople-NumOfUnemp);
    ThetaDraws = [0*ones(1,NumOfUnemp),ThetaDraws_tmp];
    PermDraws = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(1),NumOfPeople);
    for i = 1:NumOfPeriodsToSimulate
        if i==1
            ThetaList(i,:) = ones(1,NumOfPeople);
            PermList(i,:) = ones(1,NumOfPeople);
        else
            ThetaList(i,:) = randsample(ThetaDraws,NumOfPeople);
            PermList(i,:) = randsample(PermDraws,NumOfPeople);
        end
    end
else                                                                       % If st. dev. of income is age-varying, then simulation becomes slow
    for i = 1:NumOfPeriodsToSimulate
        if i<=41                                                           % Draws for pre-retirement years
            if i==1                                                        % At s=25, all draws equal 1
                ThetaDraws = ones(1,NumOfPeople);
                PermDraws = ones(1,NumOfPeople);
            else
                ThetaDraws_tmp = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfPeople-NumOfUnemp);
                ThetaDraws = [0*ones(1,NumOfUnemp),ThetaDraws_tmp];
                PermDraws = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfPeople);
            end
        else                                                               % Draws for post-retirement years
            error('Simulating retirement years as well. Write the code for retirement year draws.')
        end
        ThetaList(i,:) = randsample(ThetaDraws,NumOfPeople);
        PermList(i,:) = randsample(PermDraws,NumOfPeople);
    end
end
clear ThetaDraws_tmp ThetaDraws PermDraws