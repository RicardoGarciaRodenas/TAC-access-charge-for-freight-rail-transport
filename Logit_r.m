function [RailwayShareMarket] = Logit_r(o,d,parameters,A_r,lambda_r,tau_r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Share markets as function of the costs
% See equations (23)-(29) of the paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cost road mode 
V_omega=parameters.betaTruck*parameters.CostTruck*1.18 +...
    parameters.alfa(o)+parameters.alfa(d);
% cost railway mode
U_omega=parameters.betaRail* ((A_r+lambda_r)* tau_r+ ...
    parameters.CostRailKilometer);

% Railway Share Market using a logit model
RailwayShareMarket=exp(U_omega)/(exp(V_omega)+exp(U_omega));
end