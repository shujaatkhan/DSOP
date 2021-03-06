%
% Import all paramters from params.json
%{

    < Notes for code in progress: 
      Setup files to "capture" in this process:


        setup_everything
        setup_GList;                        % Load GList
        setup_Betacorr;                     % Load corrected beta
        NumOfPeriodsToSimulate = 90-55+1;   % Length of life in simulation (simulate until age 60)
        NumOfPeople            = 10000;      % Number of people to simulate
        NumOfPeopleBootstrap   = 10000;      % Number of people to simulate in bootstopping process
        NumOfBootstrap         = 30;        % Number of times to iterate bootstapping process
        setup_shocksLists                   % Set up shock lists


        % Estimation
        Data_SCF_wealth                          % Load data
        WealthCollege = WealthPopulationCollege; % Default WealthCollege is population itself
        weight        = ones(1,7);               % Weight = 1 for each of the 7 age groups
        x0            = [4.0,0.99];
        options=optimset('Display','final','MaxFunEvals',10000,'MaxIter',10000,'tolx',0.01,'tolfun',1,'OutputFcn', @outfun); 

      
    ----------------------------------------------------------------------------
    NOTE: to use the JSON load function, you need to download the "jsonlab" 
    toolbox from the mathworks matlabcentral website, here:

       http://www.mathworks.com/matlabcentral/fileexchange/33381-jsonlab--a-toolbox-to-encode-decode-json-files-in-matlab-octave

    Download the latest version of the jsonlab zip file.
    Unzip the file (use "unzip -a downloaded_file_name.zip" for Linux) and place 
    it in your MATLAB toolbox directory. 
    Finally, open MATLAB and add the "jsonlab" directory, inside the toolbox 
    directory, to the MATLAB path (or take equivalent steps in Octave).
%}








% Set up file path and 
params_file_path = '../../Parameters/params.json'
    % Note: this assumes that MATLAB is running in the 
    %   'SolvingMicroDSOPs/Code/Matlab/Structural Estimation/'
    % directory. 

json_params = loadjson(params_file_path)

#------------------------------------------------------------
% replicate guts of:   setup_everything

# Set up some global paramters:
global PeriodsToSolve... 
    RFree n InitialWYRatio InitialWYRatioProb...
    VarInitialLogInc ProbOfAlive LevelAdjustingParameter...
    ThetaMat ThetaMatProb PermMat PermVecProb...
    AlphaVec nP uP GList Betacorr NumOfPeriodsToSimulate NumOfPeople...
    ThetaList PermList stIndicator pi WealthCollege weight Constrained


   
setup_params

setup_functions

setup_grid_eee

setup_shocks








