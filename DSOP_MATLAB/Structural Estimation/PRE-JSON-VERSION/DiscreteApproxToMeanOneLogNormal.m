
% DiscreteApproxToMeanOneLogNormal.m
% Function which is used for constructing shock vectors (PermVec, TranVec, PermShockDraws, TranShockDraws) 
 
function shocklist= DiscreteApproxToMeanOneLogNormal(std,numofshockpoints)

global LevelAdjustingParameter sigma 
 % need to declare global variables since these are used in function FuncToIntegrate

LevelAdjustingParameter = -(1/2)*(std)^2;
for i=1:numofshockpoints-1
    ListOfEdgePoints(i) = logninv(i/numofshockpoints,LevelAdjustingParameter,std);
end

sigma = std; 

shocklist(1) = integral(@FuncToIntegrate,0,ListOfEdgePoints(1),'RelTol',1.e-12,'AbsTol',1.e-12)*numofshockpoints;
for i=2:numofshockpoints-1
    shocklist(i) = integral(@FuncToIntegrate,ListOfEdgePoints(i-1),ListOfEdgePoints(i),'RelTol',1.e-12,'AbsTol',1.e-12)*numofshockpoints;
end
shocklist(numofshockpoints) = integral(@FuncToIntegrate,ListOfEdgePoints(numofshockpoints-1),inf,'RelTol',1.e-12,'AbsTol',1e-12)*numofshockpoints;

% Note: the quad function is being phased out of Matlab in favor of
% integral, which allows the user to specify the relative and absolute
% tolerances individually, as well as an explicit use of the indefinite 
% upper limit on integration. This change of numerical intrgrator, as well 
% as specifying tight tolerlances on the allowed error, allows for a higher
% degree of accuracy in the discrete approximation than what was previosuly
% attainable under the quand function. An appropriate way to test the
% accuracy of the approximation is to examine the resulting expected value,
% discrete_shocks*discrete_probs, for which we know the true analytical 
% solution. 

%shocklist(1) = quad(@FuncToIntegrate,0,ListOfEdgePoints(1))*numofshockpoints;
%for i=2:numofshockpoints-1
%    shocklist(i) = quad(@FuncToIntegrate,ListOfEdgePoints(i-1),ListOfEdgePoints(i))*numofshockpoints;
%end
%shocklist(numofshockpoints) = quad(@FuncToIntegrate,ListOfEdgePoints(numofshockpoints-1),10)*numofshockpoints;
