%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Baseline Parameter values                                                %
%                                                                          %
%   Outputs:                                                               %
%       Rho - Coefficient of Relative Risk Aversion                        %
%       Beta - Discount factor                                             %
%       NumOfThetaShockPoints - No. of points in discrete approx. to lognormal dist
%       Sigma - Standard deviation of lognormal distribution               %
%       RFree - Interest factor                                            %
%       Gamma - Permanent income growth factor                             %
%       mMin - Minimum point in mVec                                       %
%       mMax - Maximum point in mVec                                       %
%       mHuge - mHuge is a point so that extrapolation is not needed       %
%       NumOfmPts - Number of points in mVec                               %
%       GothicAMin - Lower bound for GothicAVec                            %
%       GothicAMax - Maximum point in GothicAVec                           %
%       GothicAHuge - GothicAHuge a point so that extrapolation is not needed
%       NumOfGothicAPts - Number of points in GothicAVec                   %
%       PeriodsSolved - Number of periods back from T for which the model has been solved 
%       Constrained - Indicator variable for whether constraint is binding %
%		MC - indicates it is the multicontrol problem being considered.	   %
%                                                                          %
%__________________________________________________________________________%

Rho     = 2;                                %(* Coefficient of Relative Risk Aversion *)
Beta    = 0.96;                             %(* Discount factor *)
NumOfThetaShockPoints = 7.;                 %(* Number of points in the discrete approximation to lognormal dist *)
Sigma   = 1;                                %(* Standard deviation of lognormal distribution *)
RFree   = 1.02;                             %(* Gross interest rate *)

Gamma   = 1.00;                             % (* Permanent income growth factor *)

mMin    = 0;                                %(* Minimum point in mVec *)
mMax    = 4.;                               %(* Maximum point in mVec *)
mHuge   = 5;                                %mHuge is a point so that extrapolation is not needed
NumOfmPts = 5;                              %(* Number of points in mVec *) %Perhaps change back to 5


GothicAMin  = 0.;                           %(* Lower bound for GothicAVec *)
GothicAMax  = 4.;                           %(* Maximum point in GothicAVec *)
GothicAHuge = 9000;
NumOfGothicAPts = 5;                        %(* Number of points in GothicAVec *)

PeriodsSolved = 0;                          %(* Number of periods back from T for which the model has been solved *)

Constrained = 0;                            % Constrained

# IMPORTANT NOTE: As of 23 June 2014, the *constrained* version of the Matlab
# code has some bugs. This has not yet been resolved, and the user is encouraged 
# to be cautious in using the code with Constrained = 1.  
# However, the unconstrained code is working correctly, and the user is reminded
# that when shocks to income include a 0 value, consumers will self-impose a 
# borrowing constraint, acting as if Constrained = 0.

MC = 0;                                     % Indicates if if is the multicontrol problem



