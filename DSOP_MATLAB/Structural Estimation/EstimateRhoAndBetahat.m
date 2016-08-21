% EstimateRhoAndBetahat.m
% This file estimates params based on SCF data on wealth permanent income 
% ratio, using method of simulated moment

% ------------------- NOTE: following taken from doall.m ------------------
clc;
clear;
close all;
fclose('all');
thispath = mfilename('fullpath');
thisfile = mfilename;
pathlength = size(thispath,2);
namelength = size(thisfile,2);
CurDir = thispath(1:(pathlength-namelength));
cd(CurDir);
path(CurDir,path);
path([CurDir '..'],path);
% ------------------- NOTE: above taken from doall.m ----------------------
TimeS = cputime;
disp('=======================')
%--------------------------------------------------------------------------
% Setup
% SAK 08/10/16:
% Set ind_mom to true to use Method of Moderation (MoM).
% Default is endogenous gridpoint without method of moderation.
% If method of moderation is used, then the stardard deviation of the
% transitory and permanent shocks cannot be 0. The code does not produce
% estimates. CHECK WHY. The Mathematica code that uses MoM produces an
% error if a standard deviation of 0 is used post retirement.
%--------------------------------------------------------------------------
% Under default setting, without bootstrap, program takes about 0.47mins 
% to run on a 2.7 GHz Intel Core i5 MacBook Pro.
%--------------------------------------------------------------------------
ind_mom = false;                                                           % Use MoM (Default: false)
ind_sdtv = false;                                                          % Use time-varying SD of shocks during working age (Default: false)
if ind_mom
    disp('Method of Moderation is being used.')
else
    disp('Method of Endogenous Gridpoints is being used, without MoM.')
end
disp('=======================')
if ind_sdtv
    disp('Age-varying standard deviation of shocks during working age.')
else
    disp('Fixed standard deviation of shocks during working age.')
end
disp('=======================')
%--------------------------------------------------------------------------
rng(111)                                                                   % Setting a seed for the random number generator
setup_everything
%--------------------------------------------------------------------------
setup_GList;                                                               % Load GList
setup_Betacorr;                                                            % Load corrected beta
NumOfPeriodsToSimulate = 90-55+1;                                          % Length of life in simulation (simulate from 25-60) NOTE: If we don't have shocks for s=25, then why do we need to simulate 25?
NumOfPeople            = 1000;                                             % Number of people to simulate
NumOfPeopleBootstrap   = 1000;                                             % Number of people to simulate in bootstopping process
NumOfBootstrap         = 20;                                               % Number of times to iterate bootstapping process
TimeS = cputime;
setup_shocksLists                                                          % Set up shock lists (Shocks used in simulation)
%--------------------------------------------------------------------------
% Estimation
%--------------------------------------------------------------------------
Data_SCF_wealth                                                            % Load data
WealthCollege = WealthPopulationCollege;                                   % Default WealthCollege is population itself
weight        = ones(1,7);                                                 % Weight = 1 for each of the 7 age groups
x0            = [4.0,0.99];
%options=optimset('Display','final','MaxFunEvals',10000,'MaxIter',10000,'tolx',0.01,'tolfun',1,'OutputFcn', @outfun);
options=optimset('Display','final','MaxFunEvals',10000,'MaxIter',10000,'tolx',0.01,'tolfun',1);
disp('Starting the minimization routine...')
ParamsWithPopulation                                                       % Estimate params with population data
%--------------------------------------------------------------------------
% Default results:
% Rho,      Betahat bl
% 4.4408    1.0301
%--------------------------------------------------------------------------
%Bootstrap                                                                 % Estimate params and standard errors by bootstrapping (uncomment to run)
%--------------------------------------------------------------------------
% Display time spent
TimeE = cputime;
disp('Time Spent (mins)')
disp((TimeE - TimeS)/60)
%--------------------------------------------------------------------------
