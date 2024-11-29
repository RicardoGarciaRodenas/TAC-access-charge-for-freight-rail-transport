%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x]=NumNode(node)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch node
    case {'ES'}
        x=1;
    case {'FR'}
        x=2;
    case {'IT'}
        x=3;
    case {'SL'}
        x=4;
    case {'HR'}
        x=5;
    case {'HU'}
        x=6;
    case {'SE'}
        x=7;
    case {'NE'}
        x=8;
    case {'W'}
        x=9;
end
end