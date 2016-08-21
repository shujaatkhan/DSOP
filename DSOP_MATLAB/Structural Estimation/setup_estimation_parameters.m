
% First load all values from the JSON file. NOTE that, because dont have the libary installed, will need to guess at some of the functionality.
% Two main options:  (1) everything is directly loaded into the workspace  (2)  everything is loaded into a single object, which has cells referenced
% Ah, seems like the first is the answer.

% Load all the parameters:
param = loadjson('../../Parameters/params_json_replication.json')

% Now need to set values to match the proper names from the files. Will be
% stepping through a number of individual files to confirm that we are either
% (a) pulling parameters into the corret file, or (b) properly re-creating the
% manipulations which need to happen there -- either *here* or elsewhere in
% their own seperate files.
% NOTE that may need / want to compare some of the detailed manipulations needed
% to the Py code, and aim to make both structures match for ease of later usage.
% (When possible make *Py* conform to *Mat* style.)

% Stepping though:     setup_params.m:
Constrained       = param.constrained     % 0;  % todo Note: need to comfirm that the string "false" will convert appropriately.
PeriodsToSolve    = param.final_age - param.initial_age    % 90-25; % Number of periods to iterate
NumOfShockPoints  = param.psi_N     %8; % Number of points in the discrete approximation to lognormal dist

if param.psi_N != param.xi_N
  disp('WARNING: param.psi_N != param.xi_N when reading parameters in from JSON file.')
end

% todo [NOTE][NEED TO CHECK on why the following was done in the matlab setup_params.m code. NPalmer 13Feb2015]
% VarPerm = zeros(1,PeriodsToSolve-(param.final_age - param.retirement_age));   % zeros(1,PeriodsToSolve-(90-65))
% VarTran = zeros(1,PeriodsToSolve-(param.final_age - param.retirement_age)); %  zeros(1,PeriodsToSolve-(90-65))
%
% for i = 1:length(VarPerm)
% if i<=20
%     VarPerm(i) = 0.075-(0.075-0.03)/(PeriodsToSolve-(param.final_age - param.retirement_age)-20-1)*(i-1);
% else
%     VarPerm(i) = 0.03+(0.08-0.03)/20*(i-20);
%     end
%     VarTran(i) = 0.14-(0.14-0.03)/(PeriodsToSolve-(param.final_age - param.retirement_age)-1)*(i-1);
% end
% SIGMAPerm          = [sqrt(VarPerm) 0.001*ones(1,param.final_age - param.retirement_age)];   %x
% SIGMATran          = [sqrt(VarTran) 0.001*ones(1,param.final_age - param.retirement_age)];   %x

% todo [NOTE] Note from setup_params.m: these parameters are from Carroll (1992), a little more conservative than Carroll and Samwick (1997)
% todo [NOTE] (Note2: unclear if for values above or below in code.)
% todo Need to confirm the "order" of the SIGMAPerm  and  SIGMATran values! They *appear* to potentially be in the wrong order? Need to confirm this.

SIGMAPerm =   param.psi_sigma      % todo: order-of-time tag
SIGMATran =   param.xi_sigma       % todo: order-of-time tag

pUnemp             = param.p_unemploy;         % 0;%0.5/100; % Probability of unemployment (when unemployed inc level is zero)
pRetire            = param.p_unemploy_retire;  % pUnemp/10;   %x
RFree              = param.R;                  % 1.03;  % Gross interest rate
AlphaMin           = param.a_min;              % 0.00001; % Minimum point in AlphaVec (glid of possible saving)
AlphaMax           = param.a_max;              % 4;  % Maximum point in AlphaVec
AlphaHuge          = param.a_huge;             % 10; % Value of Alpha at which we connect to perf foresight function
n                  = param.a_size;             % 50; % Number of points in AlphaVec

InitialWYRatio     = param.initial_wealth_income_ratio_vals;  %[0.17,.5,.83]; % Initial wy ratio (from the program on the paper (p.13))
InitialWYRatioProb = param.initial_wealth_income_ratio_probs; %[.33333, .33333, .33334];    % Prob associated with initial wy ratio
VarInitialLogInc   = param.variance_initial_log_income;   % Variance of initial log income
% todo: need to confirm how these three values above are used in the estimation process, and what the best approaches are

% Probability of being alive after retirement
ProbOfAlive = param.survival_probs;   % todo: order-of-time tag
% [9.8438596e-01   9.8438596e-01   9.8438596e-01   9.8438596e-01   9.8438596e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.7567062e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.6207901e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   9.3721595e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01   6.3095734e-01]; % (1st element is the prob of being alive until age 66)
% ProbOfAlive = [ones(1,PeriodsToSolve-length(ProbOfAlive)),ProbOfAlive];

% todo: confirm, again, that this is in the right order. NOTE AGAIN that the indexing may be complex in the code itself.
% todo: simple fix, of course, if there are complex time-indices in the code: simply reverse all inputs here.

% Position of match
pi = param.quantile_to_match             % 0.50; % pith quantile is matched






# ======================== THERE BE PYTHON CODE BELOW! ======================= #
from __future__ import division
from time import time
import numpy as np
import pandas as pd
import warnings
import json
import csv



# Open the json params file and read in:
infile = open("params.json",'rb')

scf_data_path = './'

json_params = json.load(infile)
infile.close()

# ------------------------------------------------------------------------------
# -------------- Set up decision problem and simulation parameters -------------
# ------------------------------------------------------------------------------

# SET A SEED:
seed = int(json_params['seed'])
np.random.seed(seed)
    # NOTE: need to determine a better way to do this, instead of changing the
    # default/global random seed. Ideally we would be able to use a
    # numpy.random.RandomState RNG object, but if use this, will need to change
    # all the scipy.stats functions used. Return to...
    #
    # Regardless, when running "bootstrap" process, I presume only point of
    # variation is supposed to come from the bootstrap sample *only* -- yes?
    # Confirm that we *retain* the sample income shocks drawn -- essentially,
    # only add variance via the bootstrap.
    #
    # Need to think about a little more -- and finally, ask.  (SO, etc).


bootstrap_sample_size = int(json_params["bootstrap_sample_size"])
Num_agents = int(json_params['Nagents'])
Nagents_bootstrap = int(json_params['Nagents_bootstrap'])
spline_k = int(json_params['spline_k'])
TT = int(json_params['total_simuation_periods_including_T'])
    # Total time: total periods including the final period.

index_t0 = TT-1     # Index of initial period. Subtract 1 due to zero-indexing;
                    # this is initial period due to the 'T-t' timing convention.
                    # We declare this independently of the various lifetime-
                    # length variable arrays (eg Gamma) to allow for (a) the
                    # researcher to independently declare the sim length, and
                    # (b) thus double-check all vector values. It is easy to
                    # make a mistake aligning all the lifetime vectors; this
                    # adds one level of double checking.

# Decision problem parameters:
rho_start = float(json_params['rho_first_guess'])
beta_start = float(json_params['beta_first_guess'])
R = float(json_params['R'])
beta_upper_bound = float(json_params['beta_upper_bound'])
beta_lower_bound = float(json_params['beta_lower_bound'])
rho_upper_bound = float(json_params['rho_upper_bound'])
rho_lower_bound = float(json_params['rho_lower_bound'])
constrained = bool(json_params['constrained'])

Gamma  = np.array(json_params['Gamma_vector'])
timevary_discount_factors = np.array(json_params['timevary_discount_factors'])
survival_probs = np.array(json_params['survival_probs'])

# Income process parameters
psi_sigma = np.array(json_params['psi_sigma'])
psi_N = int(json_params['psi_N'])
xi_sigma = np.array(json_params['xi_sigma'])
xi_N = int(json_params['xi_N'])
final_work_index = int(json_params['final_work_index'])

p_unemploy = float(json_params['p_unemploy'])
p_unemploy_retire = float(json_params['p_unemploy_retire'])

value_unemploy = float(json_params['value_unemploy'])
value_unemploy_retire = float(json_params['value_unemploy_retire'])

# Set up the a-grid:
grid_type = json_params['grid_type']
exp_nest = int(json_params['exp_nest'])
a_min = float(json_params['a_min'])
a_max = float(json_params['a_max'])
a_size = int(json_params['a_size'])
a_extra = float(json_params['a_extra'])
a_huge = json_params['a_huge']
if a_huge is not None:
    a_huge = float(a_huge)    # TODO: we will need to figure out how to handle "None" in Matlab.


# Set up values for estimation:
initial_wealth_income_ratio_vals = json_params['initial_wealth_income_ratio_vals']
initial_wealth_income_ratio_probs = json_params['initial_wealth_income_ratio_probs']

initial_wealth_income_ratio_seed = int(json_params['initial_wealth_income_ratio_seed'])  # Set a seed for the bootstrap

quantile_to_match = float(json_params['quantile_to_match'])


# ------------------------------------------------------------------------------
# ------------------- Set up the SCF wealth data --------------------------------
# ------------------------------------------------------------------------------

#TODO: Convert the following to lists and numpy arrays!

# Version which works via csv reader:

infile = open(scf_data_path + 'SCFdata.csv', 'rb')  # Open file handle
csv_reader = csv.reader(infile)                     # Create reader object

scf_header = csv_reader.next()                      # Pull the csv header. We
                                                    # will use for selecting
                                                    # the columns.

# Pull the column index for each group from the scf_header:
scf_m_col = scf_header.index('wealth_income_ratio')
scf_ages_col = scf_header.index('age_group')
scf_weights_col = scf_header.index('weight')

# Create a dictionary to hold the scf values we will read in. Initiate empty
# lists for wealth_to_income ("m"), the age/cohort group number ("ages"), and
# the weights for each obs ("weights").
# NOTE: these must match the values in the "scf_header"
scf_data_dictionary = {"wealth_income_ratio":[],
                       "age_group":[],
                       "weight":[]}

# Now read in the data from the datafile. The CSV-reader will loop over each
# record (row) in the file, and we will save the record into the appropriate
# spot in the above dictionary.
for line in csv_reader:
    # Read the element at line index 'wealth_income_index', convert to double,
    # and append to the scf_m list.
    scf_data_dictionary['wealth_income_ratio'].append(np.float64(line[scf_m_col]))

    scf_data_dictionary['age_group'].append(np.float64(line[scf_ages_col]))

    scf_data_dictionary['weight'].append(np.float64(line[scf_weights_col]))

infile.close()  # Close the file after finish using it.

# Now create a large data matrix from our data.
# First create an empty list, then fill it with the lists we've just created,
# using the column index locations indicated by the scf_header from above. This
# will ensure that the matrix has columns in the same order as the data file we
# read in. (Just to help with intuition.)
scf_data = []  # Empty list to fill with above.
for column_name in scf_header:
    scf_data.append(scf_data_dictionary[column_name])

# Now convert this into a matrix:
scf_data = np.array(scf_data).T


# TODO: WHERE AT: Here. Need to finish out the scf-imported-with-csv,
# and then should finally be ready to do comparisons.

# Now confirm that the SCF is organized with respect to
# empirical_age_group_indices:
# Use the "argsort" argument to find the array of indices that would sort the
# "ages" column
sorted_indices = np.argsort(scf_data[:,scf_ages_col])

# Now find the unique age groups:
empirical_age_group_indices = np.int32(np.unique(scf_data[:,scf_ages_col]))

# Ensure that the ages are sorted ascending:
empirical_age_group_indices_compare =empirical_age_group_indices.copy()
empirical_age_group_indices_compare.sort()

assert np.all(empirical_age_group_indices_compare == empirical_age_group_indices), "SCF data is not sorted!"


# Items to produce:



#-----------------------
# Pull in SCF data:
scf_dataframe = pd.read_csv(scf_data_path + 'SCFdata.csv')
scf_dataframe.describe()

# Ensure sorted appropriately by age:
if np.any(np.diff(scf_dataframe['age_group']) < 0):
    # Then need to sort by the age_group column:
    scf_dataframe.sort_index(by=['age_group'], ascending=[True], inplace=True)

# Pull out values to use in objective function:
scf_m = scf_dataframe['wealth_income_ratio']
scf_weights = scf_dataframe['weight']

empirical_age_group_indices = np.unique(scf_dataframe['age_group']) #Sorted unique values
empirical_age_group_individual_ID = scf_dataframe['age_group']
scf_age_group_counts = {}
for age in empirical_age_group_indices:
    scf_age_group_counts[age] = np.sum(scf_dataframe['age_group']==age)
    # NOTE: these ages must match the keys used below!





# We need to set up the age ranges for each cohort, and then convert those into
# the appropriate simulation indices.

# Pull the mapping between the cohort number and the actual age ranges of
# the subject observed in the data. These should be the cohorts we plan to use
# in the SMM structural estimation of the model. Once these have been defined we
# will use them below, along with the starting age, to construct a list of
# Python simulation matrix indices attached to each cohort number.

map_cohort_to_actual_age_range = {}  # Fill with numpy arrays.
for cohort_index, age_group in json_params['cohort_age_groups'].iteritems():
    # iteritems allows us to iterate over both the keys and values of a dict at
    # once. As with all dictionary operations, it does not respect any
    # particular ordering of the key-value pairs in the loop; it only maintains
    # the individual key-value pairings.
    map_cohort_to_actual_age_range[cohort_index] = np.array(age_group)

starting_age = min([min(ages) for ages in map_cohort_to_actual_age_range.itervalues()])
# The above line simply finds the minimum age in the "map_cohort..." variable.
'''
# Note: the following two-line loop can be pasted into the Python command line
# to print out the lines which define cohort_age_groups in
# seen above. If you need to change the age range, especially for a very large
# range extension, moderately tweaking this loop can make the job easier:
for j, i in enumerate(range(26, 60, 5)):
    print '"' + str(j+1) + '" :', '[', ','.join([str(n) for n in range(i, i+5)]), '],'

# To make it somewhat easier to input new values for the Gamma array, in the
# nice square shape, paste this into the Python command line.
j=0
step=5
gamma_list
for i in range(0,65,5):
    jlow=i
    jhi =i+step
    #print "jlow",jlow, "jhi",jhi
    print ', '.join([str(n) for n in testr[jlow:jhi]])
'''


# Now we programatically create the proper correspondence between actual ages
# observed in the data and the indices which appear in the simulated panel of
# income and wealth due to the 'T-t' timing convention:
simulation_map_cohorts_to_age_indices = {}

for cohort, age_indices in map_cohort_to_actual_age_range.iteritems():
    simulation_map_cohorts_to_age_indices[cohort] = index_t0 - (age_indices - starting_age)
    # On RHS: first normalize age indices to start at 0 for the lowest age,
    # then use the actual list index value of the first age to "reverse" the age
    # indices to align with the 'T-t' timing convention.


# Now the variable "simulation_map_cohorts_to_age_indices" should contain the
# appropriate cohort-number-to-simulated-age-index mapping.



# ------------------------------------------------------------------------------
# ------------------- Set up initial a0 vector ---------------------------------
# ------------------------------------------------------------------------------
# TODO: Get this!

# Set up a random number generator with a specific seed to the initial draw
# is replicable:
initial_a0_RNG = np.random.RandomState(initial_wealth_income_ratio_seed)

use_a0_bootstrap = False
if use_a0_bootstrap:
    age1_index = scf_dataframe['age_group']==1
    age_1_m = np.abs(scf_dataframe['wealth_income_ratio'][age1_index])
    a0_vector = initial_a0_RNG.choice(a=age_1_m, size=Num_agents, replace=True)
else:
    # Manually draw from the provided initial a0 distribution using inverse
    # transform sampling:
    Q = np.cumsum(initial_wealth_income_ratio_probs)  # Create cumulative discrete dist
    I = Q.searchsorted(initial_a0_RNG.uniform(0, 1, size=Num_agents))  # Inverse transform sampling

    a0_vector = np.array([initial_wealth_income_ratio_vals[i] for i in I])  # Create the draws
    # Initial tests imply this is working correctly.

# A final, unused alternative -- for testing: pull from 0,1, uniformly:
# a0_vector = initial_a0_RNG.uniform(0,1,size=Num_agents)   # VERY IMPORTANT NOTE: THIS IS JUST FILLER!
# a0_vector = np.random.choice(a=np.array([0.17, 0.5, 0.83]), size=Num_agents, replace=True) initial_wealth_income_ratio_vals  # To simply use even probs.

# Check if Nagents_bootstrap == Num_agents:
if not (Num_agents == Nagents_bootstrap):
    warnings.warn("In setup_estimation_parameters.py:    Num_agents == Nagents_bootstrap.")
    # TODO: As of this writing, Python version assumes and implements Num_agents == Nagents_bootstrap for the bootstrap potion.
