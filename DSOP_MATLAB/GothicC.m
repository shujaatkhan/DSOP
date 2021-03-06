function GothicC = GothicC(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT)
%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Gothic C - C given the derivative of Gothic V                            %
%                                                                          %
%   Inputs:                                                                %
%       z - derivative of Gothic V                                         %
%       Rho - coefficient of relative risk aversion                        %
%		Beta - the pure discount factor 								   %
%		RFree - interest risk factor 									   %
%		Gamma - the growth rate of human capital 						   %
%		NumOfThetaShockPoints - the number of points in the approx. 	   %
%		ThetaVals - The points of the discrete approximation. 			   %
%   Outputs:                                                               %
%       Gothic C - C given the derivative of Gothic V                      %
%                                                                          %
%__________________________________________________________________________%

GothicC = nP(GothicVa(a,Rho,BetaLife,RLife,GammaLife,ThetaVals,ThetaVecProb,PermVals,PermVecProb,PeriodsUntilT),Rho);
