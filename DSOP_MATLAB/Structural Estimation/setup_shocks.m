% setup_shocks.m
% Construct the possible values of the shock to income

%--------------------------------------------------------------------------
% Construct transitory shock matrix.
% Each row corresponds to an age between s=26 to s=90.
% Each column corresponds to a shock point at which expectations are
% evaluated.
% If pUnemp>0, then an additional 0-valued shock is also included in the
% matrix for all working ages, i.e. s=26 to s=65.
% If pRetire>0, then even in retirement there is a chance of Theta=0.
% The associated probabilites will be defined below.
%--------------------------------------------------------------------------
ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints);
ThetaMatProb = zeros(PeriodsToSolve,NumOfShockPoints);
if pUnemp>0
    ThetaMat = zeros(PeriodsToSolve,NumOfShockPoints+1);
    ThetaMatProb = zeros(PeriodsToSolve,NumOfShockPoints+1);
end
% Note: the retirement and non-retirement ages all face the same style of
% uncertainty
for i = 1:PeriodsToSolve
    if SIGMATran(i)~=0
        ThetaVec_tmp = DiscreteApproxToMeanOneLogNormal(SIGMATran(i),NumOfShockPoints);
    else
        ThetaVec_tmp = ones(1,NumOfShockPoints);
    end
    ThetaVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);
    if pUnemp>0
        if i<= PeriodsToSolve-(90-65)                                      % During working years
            ThetaVec = [0,ThetaVec_tmp/(1-pUnemp)];                        
            ThetaVecProb = [pUnemp,ThetaVecProb*(1-pUnemp)];
        else                                                               % During retirement years
            if pRetire>0                                                   
                ThetaVec = [0,ThetaVec_tmp/(1-pRetire)];                   
                ThetaVecProb = [pRetire,ThetaVecProb*(1-pRetire)];
            else
                ThetaVec = [1,ThetaVec_tmp];                               % If pRetire==0 and pUnemp>0, then during retirment
                ThetaVecProb = [0,ThetaVecProb];                           % we fill an additional space by any number, but
            end                                                            % assign it a probability of 0.
        end
    else
        ThetaVec = ThetaVec_tmp;
    end
    ThetaMat(i,:) = ThetaVec;
    ThetaMatProb(i,:) = ThetaVecProb;
end
%--------------------------------------------------------------------------
% Construct permanent shock matrix.
%--------------------------------------------------------------------------
PermMat = zeros(PeriodsToSolve,NumOfShockPoints);
for i = 1:PeriodsToSolve
    if SIGMAPerm(i)~=0
        PermVec = DiscreteApproxToMeanOneLogNormal(SIGMAPerm(i),NumOfShockPoints);
    else
        PermVec = ones(1,NumOfShockPoints);
    end
    PermMat(i,:)= PermVec;
end
PermVecProb = (1/NumOfShockPoints)*ones(1,NumOfShockPoints);
%--------------------------------------------------------------------------
clear ThetaVec_tmp ThetaVecProb PermVec                                    % clear unnecessary variables
