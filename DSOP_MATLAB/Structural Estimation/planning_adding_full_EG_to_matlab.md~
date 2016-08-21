Planning Adding EG to Matlab code 
=================================

First step is to simply walk through the current version of the Matlab code and see where this would be implemented. 

Steps in Matlab Solution
------------------------

The main estimation occurs in this file:

    ConstructcInterpFunc

The BIG Qs:

    - Does the big loop produce the vectors of values needed to create the function
            cons(m) = c?
      Because if it does, we can "simply" swap 
    
    
    - Note that the "magic" happens in the "ScriptC" function in GothicVa. We'll need to change this to the appropriate function, and I think it should work fine. 
    >>
    >>    for i=1:length(ThetaVals)
    >>        for j = 1:length(PermVals)
    >>        mtp1 = (RLife(PeriodsUntilT+1)/(GammaLife(PeriodsUntilT+1)*PermVals(j)).*a + ThetaVals(i));
    >>        GothicVa = GothicVa+((GammaLife(PeriodsUntilT+1)*PermVals(j))^(-Rho))*...
    >>                    uP(ScriptC(mtp1,mLowerBoundLife,DeltaGothicHLife,KappaMin,PeriodsUntilT,chiIntData,Constrained),Rho)...
    >>            *PermVecProb(j)*ThetaVecProb(i);
    >>        end
    >>    end
    >>
    
    - The goal above is to get the "simple EG" version of scriptC working on vectors, appropriately. This is *nearly* everything we need here, I believe.
    - It may be the case that we need to simply use the different scriptC function when executing EG, and nothing else. 
    - Then simple "if" statements, embedded in this bit of the code, should be able to call either version and output the approriate values. 
    
    - Next steps:
        * write a version of EG 
        * figure out a way to compare the two versions
            -- output two grids
            -- evaluate difference across these two grids
        * *Then* time to compare Matlab and Python.
    
    
    
    
The BIG OBSERVATIONS:

    - as in the earlier Python code, there is a lot of stuff calculated multiple 
      times in each call of the ObjectiveFunction in the Matlab code. 
      Specifically

