% ConstructcInterpFunc_EG.m
%
%  NMPalmer notes:
%   It appears that full creation of the interpolation matrices happens here.
%   The only question is whether these matrices are used as I might imagine --
%   to construct a function cons(m) -> c.
%   If that is the case then it should be relatively simple to replace guts with
%   EG-only calculations.
%   NMPalmer Sep 2014
%
% Use the method of endogenous gridpoints. The main difference between
% this and the method of moderation is the lack of use of second
% derivative values to bound the consumption function. In code, the
% functions GothicVa_EG and ScriptC_EG are where the main differences
% are located.

RLife = ones(1,PeriodsToSolve)*RFree;
BetaLife = zeros(1,PeriodsToSolve);
GammaLife = fliplr(GList);
cLife = ones(1,PeriodsToSolve+1);
GothicAVec = AlphaVec;
%--------------------------------------------------------------------------
% Interpolation data for 'c' and 'm' for periods s=90 to s=25 (inclusive).
% IntData(:,:,end) is for s=25 and IntData(:,:,1) is for s=90.
%--------------------------------------------------------------------------
IntData = zeros(length(GothicAVec),2,PeriodsToSolve+1);
IntData(:,:,1) = [(1:length(GothicAVec))' (1:length(GothicAVec))'];        % Interpolation function for T=90: consume everything.
Rho = x(1);
for l=1:PeriodsToSolve                                                     % Solving backwards from s=89 to s=25 (inclusive)
    P            = ProbOfAlive(length(ProbOfAlive)-l+1);                   % Probability of being alive until next period
    BetaLife(l)  = x(2)*Betacorr(length(Betacorr)-l+1)*P;                  % Theta_t+1, i.e. theta shock points from s=90 to s=26
    ThetaVec     = ThetaMat(length(ThetaMatProb)-l+1,:);                   
    ThetaVecProb = ThetaMatProb(length(ThetaMatProb)-l+1,:);
    PermVec      = PermMat(length(ThetaMatProb)-l+1,:);                    % Psi_t+1, i.e. perm shock points from s=90 to s=26
    GothicAVect = GothicAVec;                                              % Not dealing with the constraint in the default NoMoM case
    
    futureMvect = IntData(:,1,l);
    futureCvect = IntData(:,2,l);
    GothicVaVect = GothicVa_EG(GothicAVect,Rho,BetaLife,RLife,GammaLife,ThetaVec,ThetaVecProb,PermVec,PermVecProb,l-1,futureMvect,futureCvect,Constrained);
    
    cVect = nP(GothicVaVect,Rho);
    mVect = GothicAVect+cVect;
    
    IntData(:,:,l+1) = [mVect' cVect'];                                    % Update the main interpolation data:
    
    %     if min(ThetaVec) > 0 && Constrained == 1;
    %         GothicAVect = GothicAVec + aLowerBoundt;
    %         GothicAVect = [0,GothicAVect];
    %         GothicAVect = sort(GothicAVect);
    %     else
    %         GothicAVect = sort([0.005 GothicAVec])+ aLowerBoundt;
    %         GothicAVect = sort(GothicAVect);
    %     end;
    
end;


