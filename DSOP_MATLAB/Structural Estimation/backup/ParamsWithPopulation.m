% ParamsWithPopulation.m
% This file estimates params with population data

% lb = [0 0.01];
% ub = [10 1];
% options = optimset('Hessian',{'lbfgs',10});
% options = optimset(options,'GradConstr','off');
% options = optimset(options,'GradObj','off');
% options = optimset(options,'MaxFunEvals', 1e6);
% options = optimset(options,'Algorithm','interior-point');
% xMin = fmincon(@ObjectFunction,x0,...
% [],[],[],[],lb,ub,[],options);

[xMin,fval]   = fminsearch(@ObjectFunction,x0,options); 

% Display results
disp('Estimated params with population data')
disp('Rho, Betahat bl')
disp(xMin)