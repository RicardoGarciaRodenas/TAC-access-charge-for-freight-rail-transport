function [A_r,lambda_r]=CostRailway(r,parameters, N,pr)

%A_r=max( (parameters.CostRailTime / N.L(r) ) * ( 1-  N.v(r) / N.CommercialSpeed ),0);
A_r=(parameters.CostRailTime / N.L(r) ) * ( 1-  N.v(r) / N.CommercialSpeed );
lambda_r=( pr / N.bar_tau(r) ) * parameters.CostRailKilometer;

end