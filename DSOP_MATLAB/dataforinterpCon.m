%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% dataforinterpCon.m                                                       %
%                                                                          %
% Creates data matrix corresponding to the final period of life. They are  %
% used to store data used by interpolating functions within the solution   %
% method. 																   %
%                                                                          %
%__________________________________________________________________________%


% Try these values:
% vInterpData
% scriptcInterpData
% gothiccInterpData
% chiIntData
% koppaIntData
% CapitalKoppaIntData

% CapitalChiIntData(:,:,PeriodsSolved+1) = [muVect CapitalChiVals CapitalChimuVals];
% LambdaIntData(:,:,PeriodsSolved+1) = [mVect LambdamVectRealst];

chiIntData(:,:,1)= [zeros(length(m)+1,1) zeros(length(m)+1,1) zeros(length(m)+1,1)];     % TODO: ditto below
koppaIntData(:,:,1) = [zeros(length(m)+1,1) zeros(length(m)+1,1) zeros(length(m)+1,1)];  % TODO: ditto below
vInterpData(:,:,1) = [zeros(length(m)+1,1) zeros(length(m)+1,1)];           % TODO: check with Jiaxiong for "why these +1?
scriptcInterpData(:,:,1) = [zeros(length(m)+2,1) zeros(length(m)+2,1)];     % TODO: ditto above
gothiccInterpData(:,:,1) = [zeros(length(m)+2,1) zeros(length(m)+2,1)];     % TODO: ditto above
CapitalKoppaIntData(:,:,1) = [zeros(length(m)+1,1) zeros(length(m)+1,1) zeros(length(m)+1,1)];  % TODO: ditto above
CapitalChiIntData(:,:,1) = [zeros(length(m)+1,1) zeros(length(m)+1,1) zeros(length(m)+1,1)];    %TODO: same as above
LambdaIntData(:,:,1) = [zeros(length(m)+1,1) zeros(length(m)+1,1)];                             %TODO: same as above
