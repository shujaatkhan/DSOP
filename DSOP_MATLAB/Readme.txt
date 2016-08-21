This directory contains all of the files necessary for reproducing the figures up to the structural estimation section 
in the lecture notes for Solution Methods for Micro Dymanic Stochastic Optimization in Matlab.

To plot all figures use 'doall.m'.

The files replicate the work done by the equivalent Mathematica notebook.

--------------------------------------------------------------------------------
Important note: As of 23 June 2014, the *constrained* version of the Matlab code 
has some bugs. This has not yet been resolved, and the user is encouraged to be 
cautious in using the code with Constrained = 1. 

However, the unconstrained code is working correctly, and the user is reminded
that when shocks to income include a 0 value, consumers will self-impose a 
borrowing constraint, acting as if Constrained = 0.

(This same warning appears in the "setup_params.m" file, here:
    SolvingMicroDSOPs/Code/Matlab/setup_params.m
 When the issue is resolved, this warning is to be removed from both locations.)
    

