%==========================================================================%
% Solution Methods for Micro Dymanic Stochastic Optimization               %
%                                                                          %
% Does all the work used for the structural estimation.                    %
% 			                                                               %
%__________________________________________________________________________%
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



EstimateRhoAndBetahat;
plotWealthPermIncRatio;
plotFigA1_ProbOfALive;
plotFigA1_ProbOfALive;
plotFigA2_Beta;
plotFigA3_G; 