

setup_Betacorr % Load data 
setup_GList    % Load data
setup_params   % Load data

age = [26:1:90];

% Print out values to a file:

delim = ',';

# Beta:
X = [Betacorr'];
filename = ['Beta.csv'];
fid = fopen(filename, 'w');
fprintf(fid, ['Beta\n'])
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.60f', 'delimiter', delim);


# Gamma:
X = [GList'];
filename = ['Gamma.csv'];
fid = fopen(filename, 'w');
fprintf(fid, ['Gamma\n'])
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.60f', 'delimiter', delim);


# ProbOfAlive:
X = [ProbOfAlive'];
filename = ['ProbOfAlive.csv'];
fid = fopen(filename, 'w');
fprintf(fid, ['ProbOfAlive\n'])
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.60f', 'delimiter', delim);


# All together:
X = [age' Betacorr' GList' ProbOfAlive'];
filename = ['Beta_GList_ProbOfAlive.csv'];
fid = fopen(filename, 'w');
fprintf(fid, ['age_' delim 'Beta' delim 'Gamma' delim 'ProbOfAlive\n'])
fclose(fid);

% Write the data file
dlmwrite(filename, X, '-append', 'precision', '%.60f', 'delimiter', delim);

