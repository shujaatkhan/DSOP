% pre_ConstructcInterpFunc_EG.m
% A function to help implement EG in the estimation process. 
% 
% Largely taken from   doall.m   and   EstimateRhoAndBetahat.m
% 
% NMPalmer Sep 2014
  
  
%================
%   clean up
%================
clc;
clear;
close all;
fclose('all');
%================
%running all problems
%================
thispath = mfilename('fullpath');
thisfile = mfilename;
pathlength = size(thispath,2);
namelength = size(thisfile,2);
CurDir = thispath(1:(pathlength-namelength));
cd(CurDir);
path(CurDir,path);
path([CurDir '..'],path);
% ------------------- NOTE: above taken from doall.m ---------------------------
  
  
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
NumOfPeople            = 10000;      % Number of people to simulate
NumOfPeopleBootstrap   = 10000;      % Number of people to simulate in bootstopping process
NumOfBootstrap         = 30;        % Number of times to iterate bootstapping process
setup_shocksLists                   % Set up shock lists


x0            = [4.0,0.99];     % Define first guess for the 

x = x0;








