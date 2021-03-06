% export_SIGMA_vals.m
%
% A simple script to run after running the setup_shocks.m script. Export the 
% SIGMAPerm  and  SIGMATran  vectors to SIGMAPerm.csv and SIGMAPerm.csv, 
% respectively. These will be used to populate the Python replication JSON file.
% 
% This script assumes that setup_shocks.m (and any necessary scripts for that 
% script) have been run already.


% Assume that setup_shocks.m script has been run, such that SIGMAPerm and 
% SIGMATran exist in active memory


% Write the solution to a file

% First create single matrix to write:
X = [SIGMAPerm'];

% choose delimiter:
delim = ',';

% Open a file and write the header:
filename = 'SIGMAPerm.csv';
fid = fopen(filename, 'w');
fprintf(fid, ['SIGMAPerm\n']);
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.20f', 'delimiter', delim);


% Now the transitory component:
X = [SIGMATran'];

% choose delimiter:
delim = ',';

% Open a file and write the header:
filename = 'SIGMATran.csv';
fid = fopen(filename, 'w');
fprintf(fid, ['SIGMATran\n']);
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.20f', 'delimiter', delim);







