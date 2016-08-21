% dofirst.m
%
% This is code to execute the very first lines in "doall.m", without actually 
% executing the rest of that file. These commands are needed to set the path 
% correctly, which allows us to independetly run a number of individual scripts
% and functions. 
%
% NMPalmer Sep 2014


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

