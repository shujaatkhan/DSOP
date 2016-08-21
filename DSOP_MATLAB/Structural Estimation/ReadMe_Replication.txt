ReadMe_Replication.txt

N. Palmer

13 Feb 2015

There are multiple steps to confirming the replication here:

1. [x] Replicate by carefully tracking down and pulling out param values in Matlab, putting into JSON, and having Python pull from JSON
    * Done. Initial replications worked well enough.

2. [ ] Get Matlab to pull in JSON
    * [x] need to carefully review the ID-ed spots for all parameters. build a "paremeter-calling tree" to map code
    * [ ] slowly work through the MATLAB "parameter-calling tree" created above. For each step in the "tree," read parameters in from JSOn and eliminate the file.
        * Note: have each Python and Matlab code file use the same JSON file, here:
            SolvingMicroDSOPs_CFPB/Code/Parameters/params_json_replication.json
    * [ ] run old matlab vs new matlab and confirm replication.

3. [ ] replicate against Mathematica as well
